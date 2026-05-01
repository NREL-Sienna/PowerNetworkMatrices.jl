# Benchmark the multi-threading benefit of `KLULinSolvePool` when building
# constraints from `VirtualMODF` rows.
#
# Mirrors the `Threads.@spawn`-per-(monitored, contingency) work item pattern
# used by PowerSimulations.jl when adding PTDF/MODF flow expressions
# (see PowerSimulations.jl, src/devices_models/devices/AC_branches.jl, the
# `add_expressions!` method for `PTDFBranchFlow`). Each work item is one
# `vmodf[monitored_arc, contingency_modification]` query, which goes through
# `with_worker` on the pool.
#
# Usage (from the package root):
#     julia --project=test --threads=auto test/performance/bench_parallel_modf.jl
# Optional first argument: system name (default
# `"matpower_ACTIVSg2000_sys"`). The script auto-resolves the case-builder
# category. Other systems probed by the perf tests:
#     Base_Eastern_Interconnect_515GW   (70k buses; PSSE, large pool win)
#     matpower_ACTIVSg2000_sys          (2k buses; default, balanced)
#     c_sys14, c_sys5                   (small, parallel overhead dominates)

using AppleAccelerate
using PowerNetworkMatrices
using PowerSystems
using PowerSystemCaseBuilder
using Logging
using Printf

const PNM = PowerNetworkMatrices
const PSB = PowerSystemCaseBuilder
const PSY = PowerSystems

configure_logging(; console_level = Logging.Error)

const SYSTEM_NAME = length(ARGS) >= 1 ? ARGS[1] : "matpower_ACTIVSg2000_sys"
const N_REPEATS = 3  # median of N runs to dampen noise

# Match the system name to the right PSB category. Order is tried in turn;
# the first category that actually contains the system wins.
const SYSTEM_CATEGORIES = [
    PSB.MatpowerTestSystems,
    PSB.PSSEParsingTestSystems,
    PSB.PSITestSystems,
]

function resolve_system_category(name::AbstractString)
    for cat in SYSTEM_CATEGORIES
        try
            return cat, PSB.build_system(cat, name)
        catch
            continue
        end
    end
    error("No PowerSystemCaseBuilder category contains system '$(name)'.")
end

"""
    build_modf_with_contingencies(sys; nworkers, n_contingencies, n_monitored)
        -> (vmodf, work)

Build a `VirtualMODF` and register `n_contingencies` single-arc outages, then
generate `n_contingencies × n_monitored` `(monitored_arc, contingency_uuid)`
work items. The product is bounded so this stays tractable on large systems
(`n_arcs²` would be billions on the 60k Eastern Interconnect).
"""
function build_modf_with_contingencies(sys::PSY.System;
    nworkers::Int,
    n_contingencies::Int = 32,
    n_monitored::Int = 256,
)
    vmodf = PNM.VirtualMODF(sys; nworkers = nworkers)
    arc_count = length(vmodf.axes[1])
    n_contingencies = min(n_contingencies, arc_count)
    n_monitored = min(n_monitored, arc_count)
    for e in 1:n_contingencies
        b_e = vmodf.arc_susceptances[e]
        ctg_uuid = Base.UUID(UInt128(e))
        ctg = PNM.ContingencySpec(
            ctg_uuid,
            PNM.NetworkModification(
                "outage_arc_$e",
                [PNM.ArcModification(e, -b_e)],
            ),
        )
        vmodf.contingency_cache[ctg_uuid] = ctg
    end
    work = Tuple{Int, Base.UUID}[]
    sizehint!(work, n_contingencies * n_monitored)
    for e in 1:n_contingencies, m in 1:n_monitored
        push!(work, (m, Base.UUID(UInt128(e))))
    end
    return vmodf, work
end

# Reset row + woodbury caches so each timed run does the same compute work.
function reset_caches!(vmodf::PNM.VirtualMODF)
    PNM.clear_caches!(vmodf)
    return nothing
end

function run_serial(vmodf, work)
    results = Vector{Vector{Float64}}(undef, length(work))
    @inbounds for k in eachindex(work)
        m, uuid = work[k]
        ctg = vmodf.contingency_cache[uuid]
        results[k] = vmodf[m, ctg.modification]
    end
    return results
end

function run_parallel(vmodf, work)
    # One task per OS thread pulling chunks; `@spawn`-per-item would burn
    # tens of ms in scheduler overhead at 8k+ items.
    results = Vector{Vector{Float64}}(undef, length(work))
    Threads.@threads :dynamic for k in eachindex(work)
        m, uuid = work[k]
        ctg = vmodf.contingency_cache[uuid]
        results[k] = vmodf[m, ctg.modification]
    end
    return results
end

function median_time(f, vmodf, work; warmup::Bool = true)
    warmup && (reset_caches!(vmodf); f(vmodf, work))
    times = Float64[]
    for _ in 1:N_REPEATS
        reset_caches!(vmodf)
        t = @elapsed f(vmodf, work)
        push!(times, t)
    end
    sort!(times)
    return times[fld(N_REPEATS + 1, 2)]
end

function results_match(a, b; atol = 1e-9)
    length(a) == length(b) || return false
    @inbounds for k in eachindex(a)
        isapprox(a[k], b[k]; atol = atol) || return false
    end
    return true
end

function main()
    nthreads = Threads.nthreads()
    if nthreads <= 1
        @warn "Threads.nthreads() == 1; relaunch with `--threads=auto` or " *
              "`--threads=N` (N > 1) to observe parallel speedup."
    end

    println("System:           $(SYSTEM_NAME)")
    println("Threads:          $(nthreads)")
    println("Repeats:          $(N_REPEATS) (median)")
    println()

    cat, sys = resolve_system_category(SYSTEM_NAME)
    println("Category:         $(cat)")
    println("Buses:            $(length(PSY.get_components(PSY.ACBus, sys)))")

    # Pool defaults to nthreads()-1 to leave one thread for the caller; the
    # benchmark wants every available thread on solver workers, so explicitly
    # use nthreads().
    print("Building VirtualMODF and registering contingencies... ")
    flush(stdout)
    build_t = @elapsed (vmodf, work) =
        build_modf_with_contingencies(sys; nworkers = nthreads)
    @printf("%.2f s\n", build_t)
    n_arcs = length(vmodf.axes[1])
    println("Arcs:             $(n_arcs)")
    println("Pool workers:     $(PNM.nworkers(vmodf))")
    println("Work items:       $(length(work))")
    println()

    serial_t = median_time(run_serial, vmodf, work)
    parallel_t = median_time(run_parallel, vmodf, work)
    speedup = serial_t / parallel_t

    # Verify correctness against the serial baseline.
    reset_caches!(vmodf)
    serial_ref = run_serial(vmodf, work)
    reset_caches!(vmodf)
    parallel_check = run_parallel(vmodf, work)
    ok = results_match(serial_ref, parallel_check)

    @printf("Serial median:    %.4f s\n", serial_t)
    @printf("Parallel median:  %.4f s\n", parallel_t)
    @printf("Speedup:          %.2fx (ideal %.2fx)\n", speedup, Float64(nthreads))
    println()
    println(
        if ok
            "Correctness:      OK (parallel == serial)"
        else
            "Correctness:      FAIL (parallel != serial)"
        end,
    )

    if nthreads > 1 && speedup < 1.1
        @warn "Parallel run is not faster than serial. With $(nthreads) " *
              "threads expect ≥ 1.5x; got $(round(speedup; digits = 2))x. " *
              "Either the system is too small to amortize @spawn overhead, " *
              "the pool is not in use, or worker threads are starving."
    end

    return ok ? 0 : 1
end

exit(main())
