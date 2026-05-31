# ---------------------------------------------------------------------------
# Precompilation workload
# ---------------------------------------------------------------------------
#
# Runs a small but representative set of PNM operations during package
# precompilation so the *native* code for the two things that dominate "time
# to first model" is baked into the precompile image instead of being JIT
# compiled on the user's first call:
#
#   1. Ybus assembly      — `Ybus(sys)` and the matrices derived from it
#                           (Incidence / BA / ABA / PTDF / LODF / Virtual*),
#                           including the VirtualMODF Woodbury contingency-update
#                           path (`compute_woodbury_factors` /
#                           `apply_woodbury_correction`).
#   2. The sparse linear   — `klu_factorize` + `solve!` / `solve_sparse!` /
#      solvers (KLU)        `tsolve!` / `numeric_refactor!` for both the
#                           Float64 and ComplexF64 instantiations.
#
# Why the solver paths are driven explicitly
# -------------------------------------------
# `Ybus(sys)` only *assembles* a `ComplexF32` sparse matrix; it never
# factorizes. So building a Ybus does **not** compile a single line of the KLU
# wrapper. The first thing in a user's session that actually factors-and-solves
# is `PTDF` / `LODF` / `VirtualPTDF`, and today that call pays the full KLU
# compile cost on its critical path. To "facilitate the recompilation of the
# linear solvers" we drive `klu_factorize` and every solve/refactor entry point
# directly here, for **both**:
#
#   * `{Float64, Int}`    — used by ABA / PTDF / LODF / VirtualPTDF, and
#   * `{ComplexF64, Int}` — used by the Ward reduction (`klu_zl_*`),
#
# so the solver surface is compiled regardless of which high-level matrix the
# user happens to build first, and regardless of whether they ever touch the
# Ward path. Because `precompile.jl` is part of the module, any edit to the
# `KLUWrapper` sources invalidates the PNM cache and re-runs this workload —
# the linear-solver native code is rebuilt automatically alongside it.
#
# Pointer hygiene
# ---------------
# Every `KLULinSolveCache` built below holds raw libklu Symbolic/Numeric
# pointers that are valid only inside *this* (precompile) process. None of
# these caches escape the workload: they are local, used, and explicitly
# `Base.finalize`d, so no live C handle is ever captured into the serialized
# image (which would be a dangling pointer / SIGSEGV in a fresh session).
#
# Opt out
# -------
# The whole workload is gated on the `precompile_workload` preference (default
# `true`); see `set_precompile_workload`. Toggling it is a tracked Preferences
# change, so flipping it forces a clean PNM (re)precompile.

using PrecompileTools: @setup_workload, @compile_workload

# A small, diagonally dominant (hence non-singular) symmetric tridiagonal
# sparse matrix with `Int` (== Int64) indices — the exact element/index types
# PNM hands to KLU. Built fresh so nothing aliases module state.
function _precompile_tridiag_f64(n::Int = 6)
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    for i in 1:n
        push!(rows, i)
        push!(cols, i)
        push!(vals, 4.0)
        if i < n
            push!(rows, i)
            push!(cols, i + 1)
            push!(vals, -1.0)
            push!(rows, i + 1)
            push!(cols, i)
            push!(vals, -1.0)
        end
    end
    return SparseArrays.sparse(rows, cols, vals, n, n)
end

# Complex counterpart for the `klu_zl_*` path. Adding `im` on the diagonal
# keeps it strictly diagonally dominant in modulus (|4+i| > 2), so the
# factorization never trips `SingularException` during precompile.
function _precompile_tridiag_cplx(n::Int = 6)
    A = SparseArrays.SparseMatrixCSC{ComplexF64, Int}(_precompile_tridiag_f64(n))
    for i in 1:n
        A[i, i] += im
    end
    return A
end

# Exercise the entire KLU solver surface for one `{Tv, Int}` instantiation so
# its native code is cached. Nothing here escapes; the cache is finalized.
function _precompile_klu(A::SparseArrays.SparseMatrixCSC{Tv, Int}) where {Tv}
    cache = klu_factorize(A)               # ctor + symbolic_factor! + numeric_refactor!
    n = size(A, 1)
    solve!(cache, ones(Tv, n))             # dense single-RHS  (klu_*_solve)
    solve!(cache, ones(Tv, n, 2))          # dense multi-RHS
    tsolve!(cache, ones(Tv, n))            # transpose solve   (klu_*_tsolve)
    solve_sparse(cache, A)                 # structured sparse RHS (allocating)
    numeric_refactor!(cache, A)            # refactor          (klu_*_refactor)
    full_refactor!(cache, A)               # symbolic_refactor! + numeric_refactor!
    Base.finalize(cache)                   # release libklu handles eagerly
    return nothing
end

# A tiny, fully-connected 3-bus AC loop plus a shunt. The loop (not a tree)
# makes ABA non-trivial and PTDF/LODF well defined; one REF / one PV / one PQ
# bus exercises the bus-type handling, and the `FixedAdmittance` drives the
# shunt-assembly path. Pure PSY/IS construction — kept in `@setup_workload`
# (not `@compile_workload`) because that code belongs to PowerSystems, which
# precompiles it itself.
function _precompile_build_system()
    sys = PSY.System(100.0)
    function _bus(number, name, bustype)
        return PSY.ACBus(;
            number = number,
            name = name,
            available = true,
            bustype = bustype,
            angle = 0.0,
            magnitude = 1.0,
            voltage_limits = (min = 0.9, max = 1.1),
            base_voltage = 138.0
        )
    end
    b1 = _bus(1, "b1", ACBusTypes.REF)
    b2 = _bus(2, "b2", ACBusTypes.PV)
    b3 = _bus(3, "b3", ACBusTypes.PQ)
    for b in (b1, b2, b3)
        PSY.add_component!(sys, b)
    end
    function _line(name, from, to)
        return PSY.Line(;
            name = name,
            available = true,
            active_power_flow = 0.0,
            reactive_power_flow = 0.0,
            arc = PSY.Arc(; from = from, to = to),
            r = 0.01,
            x = 0.1,
            b = (from = 0.001, to = 0.001),
            rating = 2.0,
            angle_limits = (min = -1.0, max = 1.0)
        )
    end
    PSY.add_component!(sys, _line("l12", b1, b2))
    PSY.add_component!(sys, _line("l23", b2, b3))
    PSY.add_component!(sys, _line("l31", b3, b1))
    # Fixed shunt admittance on the PQ bus -> shunt-assembly path in Ybus.
    PSY.add_component!(sys, PSY.FixedAdmittance("sh3", true, b3, 0.0 + 0.2im))
    return sys
end

# Drive the high-level matrix constructors that make up a typical "first
# model" build. The Float64 KLU solve paths are hit transitively here
# (PTDF -> solve_sparse!, VirtualPTDF row query -> solve!); `_precompile_klu`
# above covers the rest of the solver surface and the complex path.
function _precompile_system_matrices(sys)
    Ybus(sys)
    Ybus(sys; make_arc_admittance_matrices = true)
    IncidenceMatrix(sys)
    BA_Matrix(sys)
    ABA_Matrix(sys; factorize = true)
    PTDF(sys)
    LODF(sys)
    vptdf = VirtualPTDF(sys)
    vptdf[first(get_arc_axis(vptdf)), :]      # lazy row solve -> KLU solve!
    vlodf = VirtualLODF(sys)
    vlodf[first(get_arc_axis(vlodf)), :]      # lazy row solve
    # VirtualMODF drives the Woodbury contingency-update path
    # (compute_woodbury_factors / apply_woodbury_correction) that the PTDF/LODF
    # solves above never touch, so without this the first user contingency-row
    # query pays the full Woodbury compile cost. An N-1 outage on arc 1 of the
    # 3-bus loop keeps the network connected (no islanding short-circuit).
    vmodf = VirtualMODF(sys)
    mod = NetworkModification(
        "precompile_outage",
        [ArcModification(1, -vmodf.arc_susceptances[1])]
    )
    monitored = length(vmodf.arc_susceptances) > 1 ? 2 : 1
    vmodf[monitored, mod]                     # lazy Woodbury solve
    return nothing
end

if get_precompile_workload()
    @setup_workload begin
        A_f64 = _precompile_tridiag_f64()
        A_cplx = _precompile_tridiag_cplx()
        # System construction is guarded: a PSY API change must never break
        # PNM precompilation. The solver workload below does not depend on it.
        sys = try
            _precompile_build_system()
        catch err
            @debug "PNM precompile: system build skipped" exception=err
            nothing
        end
        @compile_workload begin
            # Silence the constructors' own @info/@warn (subnetwork search,
            # units-base) so they don't clutter the precompile log. Methods
            # compiled inside the closure are still captured by PrecompileTools.
            Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
                # Critical: compile the linear solvers (both element types)
                # first and unconditionally.
                _precompile_klu(A_f64)
                _precompile_klu(A_cplx)
                if sys !== nothing
                    try
                        _precompile_system_matrices(sys)
                    catch err
                        @debug "PNM precompile: system matrices skipped" exception=err
                    end
                end
            end
        end
    end
end
