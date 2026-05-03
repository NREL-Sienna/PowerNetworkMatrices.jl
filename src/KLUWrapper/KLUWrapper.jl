"""
    KLUWrapper

A small, allocation-aware wrapper over `libklu` (provided by `SuiteSparse_jll`)
designed for the access patterns of `PowerNetworkMatrices`:

- Cache the symbolic and numeric factorizations of an SPD/asymmetric sparse
  matrix and reuse them across many solves.
- Refactor (numeric only, or full) without re-allocating.
- Solve dense and **sparse** right-hand sides without materializing N×N
  intermediates when the RHS is structurally sparse.

This module is intentionally lighter than `KLU.jl`: it owns no Julia-side
copies of the matrix values, exposes the symbolic/numeric split directly, and
binds only the SuiteSparse_long (`klu_l_*`, `klu_zl_*`) entry points used by
the package.
"""
module KLUWrapper

import LinearAlgebra
import SparseArrays
import SparseArrays: SparseMatrixCSC, getcolptr, rowvals, nonzeros

"""
    KLU_POOL_DEBUG :: Bool

Compile-time gate for the KLU wrapper's runtime diagnostics — the
precondition snapshot in `solve!` and the per-worker collision detector
in `with_worker`. When `false` (default), both are folded out by the
compiler via `@static if` and contribute zero runtime cost. Flip to
`true` when reproducing platform-specific failures (e.g. the Windows
libklu KLU_INVALID + access violation thread).

This is a `const` rather than a `Preferences.@load_preference` flag so
that toggling it forces a precompile rebuild and the production binary
never carries the diagnostic code.
"""
const KLU_POOL_DEBUG = false

export KLULinSolveCache,
    KLULinSolvePool,
    klu_factorize,
    symbolic_factor!,
    symbolic_refactor!,
    numeric_refactor!,
    full_factor!,
    full_refactor!,
    solve!,
    tsolve!,
    solve_sparse!,
    solve_sparse,
    with_worker,
    acquire!,
    release!,
    nworkers,
    n_valid,
    reset!,
    is_factored

include("klu_jll_bindings.jl")
include("klu_cache.jl")
include("solve_dense.jl")
include("solve_sparse_rhs.jl")
include("pool.jl")

end # module
