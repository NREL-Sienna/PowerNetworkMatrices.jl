# Shared dispatch for "run a solve under whatever solver we have." Lives
# outside the Virtual{PTDF, LODF, MODF} files so all three matrices share
# the same `with_solver` seam and the same KLU/AppleAccelerate factory.
#
# Windows libklu thread-safety: the MinGW build of `SuiteSparse_jll`'s
# libklu exhibits `EXCEPTION_ACCESS_VIOLATION` and `KLU_INVALID` returns
# from `klu_l_solve` under parallel solves, even with verified-distinct
# per-worker `Numeric`/`Symbolic`/`Common` pointers. The single-cache
# adapter (`KLULinSolveCache` + `solver_lock`) sidesteps the bug. We
# concentrate the platform check in `_create_klu_solver` — that is the
# only `Sys.iswindows()` site — so the parallel-pool `with_solver` method
# stays platform-agnostic.

# How many per-worker scratch buffers to allocate alongside a solver. The
# KLU pool path has one per worker so concurrent `with_worker` calls each
# see exclusive scratch; single-cache and non-pool backends (Windows
# fallback, AppleAccelerate) serialize through `solver_lock` and use a
# single scratch buffer. Defined here alongside `with_solver` and
# `_create_klu_solver` so all three Virtual{PTDF,LODF,MODF} constructors
# pick up the same accounting.
_n_scratch(K::KLULinSolvePool) = nworkers(K)
_n_scratch(_) = 1

"""
    _default_pool_workers() -> Int

Default `nworkers` for the KLU pool used by `Virtual{PTDF, LODF, MODF}`.
Returns `max(1, Threads.nthreads() - 1)` on every platform: leaves one
logical thread for the calling task, clamps to `1` when Julia is
single-threaded. Windows safety is enforced separately by
`_create_klu_solver`, which falls back to a single-cache adapter
regardless of `nworkers` and warns if `nworkers > 1` was requested.
"""
@inline _default_pool_workers() = max(1, Threads.nthreads() - 1)

"""
    with_solver(f, K, work_ba_col, temp_data, solver_lock) -> result

Acquire exclusive access to a solver and a matched per-task scratch pair
`(work_ba_col_slot, temp_data_slot)`, then invoke
`f(solver, work_ba_col_slot, temp_data_slot)`. Dispatch:

- `K::KLULinSolvePool{Float64}`: routes through `with_worker`, picking
  `work_ba_col[idx]` and `temp_data[idx]` for the acquired worker.
  Multiple tasks proceed in parallel up to `nworkers`. Constructed only
  on platforms where libklu is parallel-safe (see `_create_klu_solver`).
- `K::KLULinSolveCache{Float64}`: serializes through `solver_lock` and
  uses scratch slot `1`. The no-pool form chosen when `nworkers == 1` or
  on Windows (libklu thread-safety workaround).
- Anything else (e.g. `AppleAccelerate.AAFactorization`): serializes
  through `solver_lock` and uses scratch slot `1`.

`work_ba_col` and `temp_data` are `Vector{Vector{Float64}}` whose length
matches the pool's `nworkers` (or `1` for the lock-serialized paths). They
stay as a vector-of-vectors rather than a `Matrix{Float64}` so that each
slot is a concrete `Vector{Float64}` — `_solve_factorization` and the
AppleAccelerate extension are typed on `Vector{Float64}`, and the two
buffers have different lengths (`n_buses` vs. `n_buses - n_ref_buses`),
so they could not share one matrix anyway.
"""
function with_solver(
    f::F,
    K::KLULinSolvePool{Float64},
    work_ba_col::Vector{Vector{Float64}},
    temp_data::Vector{Vector{Float64}},
    solver_lock::ReentrantLock,
) where {F}
    return with_worker(K) do cache, idx
        f(cache, work_ba_col[idx], temp_data[idx])
    end
end

function with_solver(
    f::F,
    K::KLULinSolveCache{Float64},
    work_ba_col::Vector{Vector{Float64}},
    temp_data::Vector{Vector{Float64}},
    solver_lock::ReentrantLock,
) where {F}
    return @lock solver_lock f(K, work_ba_col[1], temp_data[1])
end

function with_solver(
    f::F,
    K::KT,
    work_ba_col::Vector{Vector{Float64}},
    temp_data::Vector{Vector{Float64}},
    solver_lock::ReentrantLock,
) where {F, KT}
    return @lock solver_lock f(K, work_ba_col[1], temp_data[1])
end

"""
    _create_klu_solver(ABA; nworkers) -> KLULinSolveCache{Float64} or KLULinSolvePool{Float64}

Factory used by `Virtual{PTDF, LODF, MODF}` constructors. Returns:

- `KLULinSolveCache{Float64}` (single-cache, lock-serialized) when
  running on Windows, or when `nworkers == 1` on any platform.
- `KLULinSolvePool{Float64}` with `nworkers` pre-factored workers
  otherwise.

This is the **only** `Sys.iswindows()` site in the package. On Windows
with `nworkers > 1` the request is overridden to a single cache and a
one-shot `@warn` informs the caller; passing `nworkers = 1` (or
running Julia single-threaded) suppresses the warning. The default
`_default_pool_workers()` returns `max(1, Threads.nthreads() - 1)` on
every platform, so accepting the default on Windows still triggers the
warning whenever Julia runs with more than one thread.
"""
function _create_klu_solver(
    ABA::SparseArrays.SparseMatrixCSC{Float64, Int};
    nworkers::Int,
)
    @static if Sys.iswindows()
        if nworkers > 1
            @warn "KLULinSolvePool with nworkers > 1 is not safe on Windows " *
                  "(SuiteSparse_jll libklu MinGW thread-safety bug). " *
                  "Falling back to a single-cache solver. Pass nworkers=1 " *
                  "to suppress this warning." nworkers maxlog = 1
        end
        return klu_factorize(ABA)
    end
    nworkers == 1 && return klu_factorize(ABA)
    return KLULinSolvePool(ABA; nworkers = nworkers)
end
