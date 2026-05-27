# SparseMultiply bindings. These do not require a factorization cache; they
# take a `SparseMatrixCSC{Float64, Int}` directly and build the Apple-side
# `SparseMatrix_t` view at the ccall boundary. Per call this allocates two
# transient arrays (0-based `Clong[]` colptr, narrowed `Cint[]` rowval); if
# a future hot loop needs allocation-free SpMM, lift these into a dedicated
# `AASparseView` cache type.

"""
    aa_spmm!(Y, A, X) -> Y

Compute `Y ← A · X` in place, where `A` is a `SparseMatrixCSC{Float64, Int}`,
`X` and `Y` are `StridedMatrix{Float64}` of compatible shape. libSparse
overwrites `Y` (does not accumulate).
"""
function aa_spmm!(
    Y::StridedMatrix{Cdouble},
    A::SparseMatrixCSC{Float64, Int},
    X::StridedMatrix{Cdouble},
)
    size(A, 2) == size(X, 1) || throw(
        DimensionMismatch(
            "A is $(size(A)), X is $(size(X)); inner dimensions must match.",
        ),
    )
    size(Y, 1) == size(A, 1) && size(Y, 2) == size(X, 2) || throw(
        DimensionMismatch(
            "Y has size $(size(Y)); expected $((size(A, 1), size(X, 2))).",
        ),
    )
    stride(X, 1) == 1 || throw(ArgumentError("X must have unit stride in dim 1."))
    stride(Y, 1) == 1 || throw(ArgumentError("Y must have unit stride in dim 1."))
    size(X, 2) == 0 && return Y
    cs, ri = _csc_to_apple(A)
    sp = SparseMatrix_t(
        SparseMatrixStructure(
            Cint(size(A, 1)),
            Cint(size(A, 2)),
            pointer(cs),
            pointer(ri),
            ATT_ORDINARY,
            UInt8(1),
        ),
        pointer(nonzeros(A)),
    )
    GC.@preserve cs ri A X Y _sparse_multiply_matrix!(
        sp,
        _dense_matrix(X),
        _dense_matrix(Y),
    )
    return Y
end

"""
    aa_spmv!(y, A, x) -> y

Compute `y ← A · x` in place.
"""
function aa_spmv!(
    y::StridedVector{Cdouble},
    A::SparseMatrixCSC{Float64, Int},
    x::StridedVector{Cdouble},
)
    size(A, 2) == length(x) || throw(DimensionMismatch(
        "A is $(size(A)), length(x) = $(length(x)).",
    ))
    length(y) == size(A, 1) || throw(DimensionMismatch(
        "length(y) = $(length(y)); expected $(size(A, 1)).",
    ))
    stride(x, 1) == 1 || throw(ArgumentError("x must have unit stride."))
    stride(y, 1) == 1 || throw(ArgumentError("y must have unit stride."))
    cs, ri = _csc_to_apple(A)
    sp = SparseMatrix_t(
        SparseMatrixStructure(
            Cint(size(A, 1)),
            Cint(size(A, 2)),
            pointer(cs),
            pointer(ri),
            ATT_ORDINARY,
            UInt8(1),
        ),
        pointer(nonzeros(A)),
    )
    GC.@preserve cs ri A x y _sparse_multiply_vector!(
        sp,
        _dense_vector(x),
        _dense_vector(y),
    )
    return y
end

# Convert CSC's 1-based `colptr::Vector{Int}` and `rowval::Vector{Int}` into
# Apple's 0-based `Clong[]` colstarts and narrowed `Cint[]` rowindices.
function _csc_to_apple(A::SparseMatrixCSC{Float64, Int})
    cp = getcolptr(A)
    rv = rowvals(A)
    cs = Vector{Clong}(undef, length(cp))
    ri = Vector{Cint}(undef, length(rv))
    @inbounds for k in eachindex(cp)
        cs[k] = Clong(cp[k] - 1)
    end
    @inbounds for k in eachindex(rv)
        ri[k] = Cint(rv[k] - 1)
    end
    return cs, ri
end
