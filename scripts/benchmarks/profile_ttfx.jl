# profile_ttfx.jl
#
# Purpose:
#   Measure PowerNetworkMatrices "time to first result" (TTFX) in a single
#   fresh Julia process: the package load time plus the *first* call to each of
#   the matrix constructors that dominate a new session — Ybus assembly and the
#   sparse-linear-solver-backed PTDF / LODF / VirtualPTDF / VirtualLODF builds.
#
#   The first call to each constructor is compile-dominated, so this is exactly
#   what the precompile workload (`src/precompile.jl`) is meant to shrink. A
#   warm (second) call is also timed to expose the steady-state runtime cost
#   that precompilation does NOT change.
#
# Usage (from repo root, with the package project active):
#   julia --project=. scripts/benchmarks/profile_ttfx.jl [LABEL] [N_BUSES]
#
#   LABEL    tag written to ttfx_results.csv (default: "run")
#   N_BUSES  size of the synthetic ring+chord grid (default: 60)
#
# A/B recipe for evaluating the workload (holds everything else constant):
#   julia --project=. -e 'using PowerNetworkMatrices: set_precompile_workload; set_precompile_workload(false)'
#   julia --project=. scripts/benchmarks/profile_ttfx.jl workload_off
#   julia --project=. -e 'using PowerNetworkMatrices: set_precompile_workload; set_precompile_workload(true)'
#   julia --project=. scripts/benchmarks/profile_ttfx.jl workload_on
#
# Flipping `set_precompile_workload` marks PNM's cache stale, so the next
# `using PowerNetworkMatrices` would *recompile* PNM. That one-time
# precompilation is NOT session TTFX (it is amortized across all future
# sessions), so this script forces it via `Pkg.precompile()` BEFORE — and
# outside of — the timed `using` block below. That keeps `T_LOAD` a pure
# image-load measurement and makes a single run per label directly comparable;
# no throwaway warm-up run is needed.

const LABEL = length(ARGS) >= 1 ? ARGS[1] : "run"
const N_BUSES = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 60

# ── Separate precompilation from the timed load ─────────────────────────────
# `Pkg.precompile` rebuilds a stale cache in worker subprocesses WITHOUT
# loading the package into this process, so the compile cost lands here
# (untimed) instead of inside the `using` below. After a `set_precompile_*`
# flip this is where the workload actually recompiles.
import Pkg
println("Ensuring PowerNetworkMatrices is precompiled (untimed) …")
Pkg.precompile("PowerNetworkMatrices")

# ── Package load time (image load only; cache is already warm) ───────────────
const _T0 = time_ns()
using PowerNetworkMatrices
import PowerSystems as PSY
const T_LOAD = (time_ns() - _T0) / 1e9

# Quiet the units-base / island / subnetwork logs (Base.CoreLogging avoids a
# Logging-stdlib dependency in the active project).
Base.CoreLogging.disable_logging(Base.CoreLogging.Warn)

const PNM = PowerNetworkMatrices

# ── Synthetic system builder (no PowerSystemCaseBuilder dependency) ──────────
# A ring of `n` buses with a few chords so the network is meshed (PTDF/LODF are
# well defined and ABA is non-trivial). Bus 1 is the single reference.
function build_grid(n::Int)
    sys = PSY.System(100.0)
    buses = PSY.ACBus[]
    for i in 1:n
        bt = i == 1 ? PSY.ACBusTypes.REF :
             (isodd(i) ? PSY.ACBusTypes.PV : PSY.ACBusTypes.PQ)
        b = PSY.ACBus(;
            number = i,
            name = "b$i",
            available = true,
            bustype = bt,
            angle = 0.0,
            magnitude = 1.0,
            voltage_limits = (min = 0.9, max = 1.1),
            base_voltage = 138.0
        )
        PSY.add_component!(sys, b)
        push!(buses, b)
    end
    function mkline(name, f, t)
        PSY.Line(;
            name = name,
            available = true,
            active_power_flow = 0.0,
            reactive_power_flow = 0.0,
            arc = PSY.Arc(; from = f, to = t),
            r = 0.01,
            x = 0.1,
            b = (from = 0.001, to = 0.001),
            rating = 5.0,
            angle_limits = (min = -1.0, max = 1.0)
        )
    end
    for i in 1:n
        j = i == n ? 1 : i + 1
        PSY.add_component!(sys, mkline("ring_$(i)_$(j)", buses[i], buses[j]))
    end
    # A handful of chords across the ring to add mesh structure.
    for k in 1:max(1, n ÷ 10)
        i = k
        j = mod1(i + n ÷ 2, n)
        i == j && continue
        PSY.add_component!(sys, mkline("chord_$(i)_$(j)", buses[i], buses[j]))
    end
    return sys
end

timed(f) = (GC.gc(); local t = time_ns(); f(); (time_ns() - t) / 1e9)

println("Building synthetic system ($(N_BUSES) buses) …")
const SYS = build_grid(N_BUSES)
const NB = length(collect(PSY.get_components(PSY.ACBus, SYS)))

# ── First-call (cold) and second-call (warm) timings ────────────────────────
results = Tuple{String, Float64}[]
push!(results, ("load (using PNM)", T_LOAD))

# Ybus — the headline path.
push!(results, ("Ybus #1 (cold)", timed(() -> Ybus(SYS))))
push!(results, ("Ybus #2 (warm)", timed(() -> Ybus(SYS))))
push!(
    results,
    ("Ybus arc-adm #1", timed(() -> Ybus(SYS; make_arc_admittance_matrices = true)))
)

# Supporting matrices.
push!(results, ("IncidenceMatrix #1", timed(() -> IncidenceMatrix(SYS))))
push!(results, ("BA_Matrix #1", timed(() -> BA_Matrix(SYS))))
push!(results, ("ABA_Matrix(factorize) #1", timed(() -> ABA_Matrix(SYS; factorize = true))))

# Linear-solver-backed: PTDF (solve_sparse!), LODF.
push!(results, ("PTDF #1 (cold)", timed(() -> PTDF(SYS))))
push!(results, ("PTDF #2 (warm)", timed(() -> PTDF(SYS))))
push!(results, ("LODF #1 (cold)", timed(() -> LODF(SYS))))

# Virtual builds + first lazy row query (KLU solve!).
let
    vptdf = VirtualPTDF(SYS)
    t_build = timed(() -> VirtualPTDF(SYS))
    push!(results, ("VirtualPTDF build #1", t_build))
    arc = first(PNM.get_arc_axis(vptdf))
    push!(results, ("VirtualPTDF row #1 (cold)", timed(() -> vptdf[arc, :])))
end
let
    vlodf = VirtualLODF(SYS)
    t_build = timed(() -> VirtualLODF(SYS))
    push!(results, ("VirtualLODF build #1", t_build))
    arc = first(PNM.get_arc_axis(vlodf))
    push!(results, ("VirtualLODF row #1 (cold)", timed(() -> vlodf[arc, :])))
end

# VirtualMODF: build + first lazy contingency-row query. The first query is the
# MODF-specific compile cost — it drives the Woodbury solve path
# (compute_woodbury_factors / apply_woodbury_correction) that nothing else above
# touches. A single-line (N-1) outage on a ring edge keeps the network meshed,
# so no islanding short-circuits the Woodbury path.
let
    vmodf = VirtualMODF(SYS)
    t_build = timed(() -> VirtualMODF(SYS))
    push!(results, ("VirtualMODF build #1", t_build))
    outaged = 1                                  # ring edge -> no islanding
    b_out = vmodf.arc_susceptances[outaged]
    mod = NetworkModification("ttfx_outage", [ArcModification(outaged, -b_out)])
    monitored = length(vmodf.arc_susceptances) > 1 ? 2 : 1
    push!(results, ("VirtualMODF row #1 (cold)", timed(() -> vmodf[monitored, mod])))
    push!(results, ("VirtualMODF row #2 (warm)", timed(() -> vmodf[monitored, mod])))
end

# ── Report ──────────────────────────────────────────────────────────────────
# "First model" = package load + every first/cold call (everything tagged
# "#1"); warm "#2" re-runs are excluded since they are steady-state, not TTFX.
const T_FIRST_MODEL = sum(t
for (name, t) in results if name == "load (using PNM)" ||
                            occursin("#1", name))

println("\n" * "="^64)
println("TTFX profile   label=$(LABEL)   buses=$(NB)   julia=$(VERSION)")
println("workload pref = ", PowerNetworkMatrices.get_precompile_workload())
println("="^64)
for (name, t) in results
    println(rpad(name, 30), " ", lpad(string(round(t * 1000; digits = 1)), 10), " ms")
end
println("-"^64)
println(rpad("Σ cold first-model (incl. load)", 30), " ",
    lpad(string(round(T_FIRST_MODEL * 1000; digits = 1)), 10), " ms")

open("ttfx_results.csv", "a") do io
    for (name, t) in results
        println(io, "$(LABEL),$(name),$(round(t * 1000; digits = 3))")
    end
    println(io, "$(LABEL),SUM_cold_first_model,$(round(T_FIRST_MODEL * 1000; digits = 3))")
end
println("\nAppended to ttfx_results.csv")
