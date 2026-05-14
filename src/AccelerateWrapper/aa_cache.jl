"""
A cached libSparse linear solver for repeated solves against the same
symmetric sparse matrix structure. `numeric_refactor!` and `solve!` allocate
nothing once the cache is built.

Float64 only. The factor type defaults to `SparseFactorizationLDLT` (pivoted
Bunch-Kaufman), which matches what the old `AppleAccelerateExt.jl` used.

`reuse_symbolic` controls whether `symbolic_refactor!` keeps the analysis;
`check_pattern` adds a structural-equality check on refactor calls and is
only consulted when reusing.
"""
mutable struct AAFactorCache
    # Apple-side 0-based, narrower-integer copies of the input CSC pattern.
    # Strip to lower triangle once during `symbolic_factor!`; subsequent
    # `numeric_refactor!` calls reuse these arrays as-is.
    columnStarts::Vector{Clong}
    rowIndices::Vector{Cint}
    nzval::Vector{Cdouble}
    n::Int
    nnz_tri::Int
    factorization_type::SparseFactorization_t
    symbolic::SparseOpaqueSymbolicFactorization
    numeric::SparseOpaqueFactorization_t
    reuse_symbolic::Bool
    check_pattern::Bool
    # Bounded reusable scratch for `solve_sparse!`. Lazy-grown on first call.
    scratch::Matrix{Cdouble}
    col_map::Vector{Int}
end

@inline _dim(cache::AAFactorCache) = cache.n

Base.size(cache::AAFactorCache) = (cache.n, cache.n)
Base.size(cache::AAFactorCache, d::Integer) = d <= 2 ? cache.n : 1
Base.eltype(::Type{AAFactorCache}) = Cdouble

"""
    is_factored(cache::AAFactorCache) -> Bool

`true` when `cache` holds a valid numeric factorization ready for `solve!` /
`solve_sparse!`. `false` after construction (before `full_factor!`) or after
the libSparse handles have been finalized.
"""
function is_factored(cache::AAFactorCache)
    return cache.numeric.status == SparseStatusOk && cache.symbolic.status == SparseStatusOk
end

"""
    AAFactorCache(A; reuse_symbolic=true, check_pattern=true,
                  factorization_type=SparseFactorizationLDLT)

Build a cache for the symmetric sparse matrix `A`. Allocates the Apple-side
structural arrays (`columnStarts`, `rowIndices`, `nzval`) sized to the lower
triangle of `A` but does **not** factorize. Call `full_factor!` (or
`symbolic_factor!` followed by `numeric_refactor!`) before `solve!`.

A finalizer frees libSparse handles on GC; call `Base.finalize(cache)` to
release them eagerly.
"""
function AAFactorCache(
    A::SparseMatrixCSC{Float64, Int};
    reuse_symbolic::Bool = true,
    check_pattern::Bool = true,
    factorization_type::SparseFactorization_t = SparseFactorizationLDLT,
)
    n = size(A, 1)
    n == size(A, 2) || throw(DimensionMismatch("matrix must be square; got $(size(A))"))

    nnz_tri = _count_lower_triangle(A)
    cache = AAFactorCache(
        Vector{Clong}(undef, n + 1),
        Vector{Cint}(undef, nnz_tri),
        Vector{Cdouble}(undef, 0),
        n,
        nnz_tri,
        factorization_type,
        _null_symbolic(),
        _null_factorization(),
        reuse_symbolic,
        check_pattern,
        Matrix{Cdouble}(undef, 0, 0),
        Int[],
    )
    _populate_lower_triangle_pattern!(cache, A)
    finalizer(_free_handles!, cache)
    return cache
end

# Count nonzeros in the lower triangle (i >= j) of A.
function _count_lower_triangle(A::SparseMatrixCSC{Float64, Int})
    cnt = 0
    rowval = rowvals(A)
    @inbounds for j in 1:size(A, 2)
        for p in nzrange(A, j)
            rowval[p] >= j && (cnt += 1)
        end
    end
    return cnt
end

# Strip A to its lower triangle and fill `cache.columnStarts` / `cache.rowIndices`
# with 0-based indices. Caller must have sized these to (n+1) and nnz_tri.
function _populate_lower_triangle_pattern!(
    cache::AAFactorCache,
    A::SparseMatrixCSC{Float64, Int},
)
    rowval = rowvals(A)
    n = cache.n
    pos = 0
    cache.columnStarts[1] = 0
    @inbounds for j in 1:n
        for p in nzrange(A, j)
            if rowval[p] >= j
                pos += 1
                cache.rowIndices[pos] = Cint(rowval[p] - 1)
            end
        end
        cache.columnStarts[j + 1] = Clong(pos)
    end
    return cache
end

# Snapshot the lower-triangle values into `cache.nzval`, growing if needed.
function _populate_lower_triangle_values!(
    cache::AAFactorCache,
    A::SparseMatrixCSC{Float64, Int},
)
    if length(cache.nzval) != cache.nnz_tri
        resize!(cache.nzval, cache.nnz_tri)
    end
    rowval = rowvals(A)
    nzv = nonzeros(A)
    pos = 0
    @inbounds for j in 1:size(A, 2)
        for p in nzrange(A, j)
            if rowval[p] >= j
                pos += 1
                cache.nzval[pos] = nzv[p]
            end
        end
    end
    return cache
end

# Pattern-match guard for `numeric_refactor!`: assert the incoming CSC's
# lower-triangle pattern is identical to what was analyzed. Same shape as
# KLU's `_check_pattern_match`, just without the off-by-one increment dance
# (Apple's arrays are 0-based but we always store them 0-based — no flipping).
function _check_pattern_match(
    cache::AAFactorCache,
    A::SparseMatrixCSC{Float64, Int},
    op::AbstractString,
)
    n = cache.n
    if size(A, 1) != n || size(A, 2) != n
        throw(DimensionMismatch("Cannot $op: cache is $(n)×$(n) but A is $(size(A))."))
    end
    rowval = rowvals(A)
    pos = 0
    @inbounds for j in 1:n
        cache.columnStarts[j] == Clong(pos) || return _pattern_mismatch(op)
        for p in nzrange(A, j)
            if rowval[p] >= j
                pos += 1
                pos > cache.nnz_tri && return _pattern_mismatch(op)
                cache.rowIndices[pos] == Cint(rowval[p] - 1) || return _pattern_mismatch(op)
            end
        end
        cache.columnStarts[j + 1] == Clong(pos) || return _pattern_mismatch(op)
    end
    pos == cache.nnz_tri || return _pattern_mismatch(op)
    return nothing
end

_pattern_mismatch(op::AbstractString) =
    throw(ArgumentError("Cannot $op: matrix has different sparsity structure."))

"""
Release the libSparse numeric and symbolic handles held by `cache`, leaving
Julia-side state intact. Idempotent.
"""
function _free_handles!(cache::AAFactorCache)
    if cache.numeric.status == SparseStatusOk
        _sparse_cleanup_factor!(cache.numeric)
        cache.numeric = _null_factorization()
    end
    if cache.symbolic.status == SparseStatusOk
        _sparse_cleanup_symbolic!(cache.symbolic)
        cache.symbolic = _null_symbolic()
    end
    return nothing
end

Base.finalize(cache::AAFactorCache) = _free_handles!(cache)

# Build the libSparse `SparseMatrixStructure` view that points into the
# cache's owned arrays. Marked symmetric + lower triangle so the
# `SparseFactorizationLDLT` path treats `cache.nzval` as the lower-triangle
# entries of a symmetric matrix.
function _structure_view(cache::AAFactorCache)
    return SparseMatrixStructure(
        Cint(cache.n),
        Cint(cache.n),
        pointer(cache.columnStarts),
        pointer(cache.rowIndices),
        ATT_SYMMETRIC | ATT_LOWER_TRIANGLE,
        UInt8(1),
    )
end

function _matrix_view(cache::AAFactorCache)
    return SparseMatrix_t(_structure_view(cache), pointer(cache.nzval))
end

"""
    symbolic_factor!(cache, A)

Free any cached symbolic/numeric factor, replace the structural arrays with
`A`'s lower-triangle pattern, and analyze. Subsequent `numeric_refactor!`
calls reuse the analysis.
"""
function symbolic_factor!(cache::AAFactorCache, A::SparseMatrixCSC{Float64, Int})
    n = cache.n
    if size(A, 1) != n || size(A, 2) != n
        throw(DimensionMismatch("Cannot factor: cache is $(n)×$(n) but A is $(size(A))."))
    end
    _free_handles!(cache)
    new_nnz_tri = _count_lower_triangle(A)
    if new_nnz_tri != cache.nnz_tri
        resize!(cache.rowIndices, new_nnz_tri)
        cache.nnz_tri = new_nnz_tri
    end
    if length(cache.columnStarts) != n + 1
        resize!(cache.columnStarts, n + 1)
    end
    _populate_lower_triangle_pattern!(cache, A)
    sym = _sparse_symbolic_factor(
        cache.factorization_type,
        _structure_view(cache),
        SparseSymbolicFactorOptions(),
    )
    sym.status == SparseStatusOk || _libsparse_throw(sym.status, "symbolic factor")
    cache.symbolic = sym
    return cache
end

"""
    numeric_refactor!(cache, A)

Refresh the numeric factor on top of the existing symbolic analysis. Errors
if `symbolic_factor!` has not been called yet.
"""
function numeric_refactor!(cache::AAFactorCache, A::SparseMatrixCSC{Float64, Int})
    cache.symbolic.status == SparseStatusOk ||
        error("AAFactorCache: call symbolic_factor! before numeric_refactor!.")
    cache.check_pattern && _check_pattern_match(cache, A, "numeric_refactor")
    _populate_lower_triangle_values!(cache, A)
    if cache.numeric.status == SparseStatusOk
        _sparse_cleanup_factor!(cache.numeric)
        cache.numeric = _null_factorization()
    end
    num = _sparse_numeric_factor(
        cache.symbolic,
        _matrix_view(cache),
        SparseNumericFactorOptions(),
    )
    num.status == SparseStatusOk || _libsparse_throw(num.status, "numeric factor")
    cache.numeric = num
    return cache
end

"""
    symbolic_refactor!(cache, A)

If `cache.reuse_symbolic`, optionally verify the structure matches and reuse
the existing analysis. Otherwise, rerun `symbolic_factor!`.
"""
function symbolic_refactor!(cache::AAFactorCache, A::SparseMatrixCSC{Float64, Int})
    if !cache.reuse_symbolic
        return symbolic_factor!(cache, A)
    end
    cache.check_pattern && _check_pattern_match(cache, A, "symbolic_refactor")
    return cache
end

"""
    full_factor!(cache, A) -> cache

Run a fresh symbolic analysis followed by a numeric factorization on `A`.
"""
function full_factor!(cache::AAFactorCache, A::SparseMatrixCSC{Float64, Int})
    symbolic_factor!(cache, A)
    numeric_refactor!(cache, A)
    return cache
end

"""
    full_refactor!(cache, A) -> cache

Refresh both factorizations on `A`. Defers to `symbolic_refactor!` (which
reuses the existing analysis when `cache.reuse_symbolic` is set) followed by
`numeric_refactor!`.
"""
function full_refactor!(cache::AAFactorCache, A::SparseMatrixCSC{Float64, Int})
    symbolic_refactor!(cache, A)
    numeric_refactor!(cache, A)
    return cache
end

"""
    aa_factorize(A; reuse_symbolic=true, check_pattern=true,
                  factorization_type=SparseFactorizationLDLT) -> AAFactorCache

Build a cache for `A` and immediately compute the full factorization.
"""
function aa_factorize(
    A::SparseMatrixCSC{Float64, Int};
    reuse_symbolic::Bool = true,
    check_pattern::Bool = true,
    factorization_type::SparseFactorization_t = SparseFactorizationLDLT,
)
    cache = AAFactorCache(
        A;
        reuse_symbolic = reuse_symbolic,
        check_pattern = check_pattern,
        factorization_type = factorization_type,
    )
    return full_factor!(cache, A)
end

"""
    _ensure_scratch!(cache, block) -> Nothing

Ensure `cache.scratch` is at least `n × block` and `cache.col_map` length
`block`. Used by `solve_sparse!`.
"""
@inline function _ensure_scratch!(cache::AAFactorCache, block::Int)
    n = cache.n
    s = cache.scratch
    if size(s, 1) != n || size(s, 2) < block
        cache.scratch = Matrix{Cdouble}(undef, n, block)
    end
    if length(cache.col_map) < block
        resize!(cache.col_map, block)
    end
    return nothing
end
