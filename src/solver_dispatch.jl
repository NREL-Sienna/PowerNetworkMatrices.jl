# Shared dispatch for "run a solve under whatever solver we have." Lives
# outside the Virtual{PTDF, LODF, MODF} files so all three matrices share
# the same `with_solver` seam and the same KLU/AppleAccelerate factory.
#
# All libklu activity in this process serializes through `_LIBKLU_LOCK`
# (see `KLUWrapper.jl`). A pool-of-independent-caches design was tried
# and removed: empirically, distinct `Numeric`/`Symbolic`/`Common` per
# thread does not prevent libklu state corruption, and the per-cache
# `solver_lock` we hold here further serializes any non-libklu work in
# the callback. One factor + one cache per Virtual matrix is what
# remains — same throughput as the pool variant once the global lock
# is in place.

"""
    with_solver(f, K, work_ba_col, temp_data, solver_lock) -> result

Acquire `solver_lock`, then invoke `f(K, work_ba_col[1], temp_data[1])`.
Two overloads: one specialized on `KLULinSolveCache{Float64}` (the KLU
backend), one generic for any other factorization (`AppleAccelerate.AAFactorization`,
etc.). Both serialize through `solver_lock`; the per-cache scratch
slot at index 1 is the only slot — `work_ba_col` and `temp_data` are
single-element vectors, kept as `Vector{Vector{Float64}}` because
`_solve_factorization` and the AppleAccelerate extension are typed on
`Vector{Float64}` and the two buffers have different lengths
(`n_buses` vs. `n_buses - n_ref_buses`).
"""
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
    _create_klu_solver(ABA) -> KLULinSolveCache{Float64}

Factory used by `Virtual{PTDF, LODF, MODF}` constructors. Always returns
a single `KLULinSolveCache`; concurrency is handled by the per-cache
`solver_lock` at the wrapper level and by `_LIBKLU_LOCK` at the libklu
boundary.
"""
function _create_klu_solver(ABA::SparseArrays.SparseMatrixCSC{Float64, Int})
    return klu_factorize(ABA)
end
