# Direct bindings into Apple's libSparse.dylib (part of the Accelerate
# framework). Only the entry points actually consumed by PowerNetworkMatrices
# are wrapped — Float64-only, no Float32 mangled aliases, no QR / Cholesky-AtA
# variants. Mangled names match what AppleAccelerate.jl uses; see
# `/Library/Developer/CommandLineTools/SDKs/MacOSX*.sdk/System/Library/
# Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/
# Versions/A/Headers/Sparse/Solve.h` for the C declarations.

const LIBSPARSE =
    "/System/Library/Frameworks/Accelerate.framework/Versions/A/" *
    "Frameworks/vecLib.framework/libSparse.dylib"

@enum SparseFactorization_t::UInt8 begin
    SparseFactorizationCholesky = 0
    SparseFactorizationLDLT = 1
    SparseFactorizationLDLTUnpivoted = 2
    SparseFactorizationLDLTSBK = 3
    SparseFactorizationLDLTTPP = 4
end

@enum SparseOrder_t::UInt8 begin
    SparseOrderDefault = 0
    SparseOrderUser = 1
    SparseOrderAMD = 2
    SparseOrderMetis = 3
    SparseOrderCOLAMD = 4
end

@enum SparseScaling_t::UInt8 begin
    SparseScalingDefault = 0
    SparseScalingUser = 1
    SparseScalingEquilibriationInf = 2
end

@enum SparseStatus_t::Int32 begin
    SparseStatusOk = 0
    SparseStatusFailed = -1
    SparseMatrixIsSingular = -2
    SparseInternalError = -3
    SparseParameterError = -4
    SparseStatusReleased = -2147483647
end

@enum SparseControl_t::UInt32 begin
    SparseDefaultControl = 0
end

# `SparseAttributes_t` is a packed bitfield in C; Julia can't express that
# cleanly, so we model it as `Cuint` and assemble the bits ourselves.
const att_type = Cuint
const ATT_TRANSPOSE = att_type(1)
const ATT_UPPER_TRIANGLE = att_type(0)
const ATT_LOWER_TRIANGLE = att_type(2)
const ATT_ORDINARY = att_type(0)
const ATT_TRIANGULAR = att_type(4)
const ATT_UNIT_TRIANGULAR = att_type(8)
const ATT_SYMMETRIC = att_type(12)

struct SparseMatrixStructure
    rowCount::Cint
    columnCount::Cint
    columnStarts::Ptr{Clong}
    rowIndices::Ptr{Cint}
    attributes::att_type
    blockSize::UInt8
end

struct SparseNumericFactorOptions
    control::SparseControl_t
    scalingMethod::SparseScaling_t
    scaling::Ptr{Cvoid}
    pivotTolerance::Float64
    zeroTolerance::Float64
end

# Defaults match Apple's SolveImplementation.h.
SparseNumericFactorOptions() = SparseNumericFactorOptions(
    SparseDefaultControl,
    SparseScalingDefault,
    C_NULL,
    0.01,
    eps(Cdouble) * 1e-4,
)

struct SparseSymbolicFactorOptions
    control::SparseControl_t
    orderMethod::SparseOrder_t
    order::Ptr{Cvoid}
    ignoreRowsAndColumns::Ptr{Cvoid}
    malloc::Ptr{Cvoid}
    free::Ptr{Cvoid}
    reportError::Ptr{Cvoid}
end

# `reportError` is fired by libSparse before it returns a failure status. The
# inline `error(unsafe_string(text))` matches AppleAccelerate.jl's pattern: it
# propagates the libSparse message as a Julia exception that unwinds back
# through the originating ccall. Avoids `@error`'s allocator/logger overhead,
# which is fragile when libSparse invokes the callback from its own threads.
# Passing libc malloc/free explicitly (the C_NULL "use Apple defaults" path is
# documented but unreliable from a non-Obj-C caller).
function SparseSymbolicFactorOptions()
    return SparseSymbolicFactorOptions(
        SparseDefaultControl,
        SparseOrderDefault,
        C_NULL,
        C_NULL,
        @cfunction(Libc.malloc, Ptr{Cvoid}, (Csize_t,)),
        @cfunction(Libc.free, Cvoid, (Ptr{Cvoid},)),
        @cfunction(text -> error(unsafe_string(text)), Cvoid, (Cstring,)),
    )
end

struct DenseVector_t
    count::Cint
    data::Ptr{Cdouble}
end

struct DenseMatrix_t
    rowCount::Cint
    columnCount::Cint
    columnStride::Cint
    attributes::att_type
    data::Ptr{Cdouble}
end

struct SparseMatrix_t
    structure::SparseMatrixStructure
    data::Ptr{Cdouble}
end

struct SparseOpaqueSymbolicFactorization
    status::SparseStatus_t
    rowCount::Cint
    columnCount::Cint
    attributes::att_type
    blockSize::UInt8
    type::SparseFactorization_t
    factorization::Ptr{Cvoid}
    workspaceSize_Float::Csize_t
    workspaceSize_Double::Csize_t
    factorSize_Float::Csize_t
    factorSize_Double::Csize_t
end

# Placeholder used to initialize `AAFactorCache.symbolic` before the first
# factor. Marked released so any accidental cleanup is a no-op.
function _null_symbolic()
    return SparseOpaqueSymbolicFactorization(
        SparseStatusReleased,
        0,
        0,
        ATT_ORDINARY,
        0,
        SparseFactorizationLDLT,
        C_NULL,
        0,
        0,
        0,
        0,
    )
end

struct SparseOpaqueFactorization_t
    status::SparseStatus_t
    attributes::att_type
    symbolicFactorization::SparseOpaqueSymbolicFactorization
    userFactorStorage::Bool
    numericFactorization::Ptr{Cvoid}
    solveWorkspaceRequiredStatic::Csize_t
    solveWorkspaceRequiredPerRHS::Csize_t
end

function _null_factorization()
    return SparseOpaqueFactorization_t(
        SparseStatusReleased,
        ATT_ORDINARY,
        _null_symbolic(),
        false,
        C_NULL,
        0,
        0,
    )
end

# Build the Apple-side dense views at the ccall boundary. `StridedMatrix`'s
# first-dimension stride must be 1 (the libSparse contract); we assert at the
# call site, not here.
function _dense_matrix(B::StridedMatrix{Cdouble})
    return DenseMatrix_t(
        Cint(size(B, 1)),
        Cint(size(B, 2)),
        Cint(stride(B, 2)),
        ATT_ORDINARY,
        pointer(B),
    )
end

function _dense_vector(b::StridedVector{Cdouble})
    return DenseVector_t(Cint(length(b)), pointer(b))
end

# --- ccalls -----------------------------------------------------------------
#
# Mangled symbol names come from the C++ ABI of libSparse. They are stable on
# the system framework and match what AppleAccelerate.jl binds. If Apple
# breaks them in a future macOS release, the failure will be loud (dlopen of
# a missing symbol at first call), which is the behavior we want.

# Symbolic-only factor: analyzes the pattern, returns an opaque symbolic
# factor. Can back many numeric factors on the same pattern.
function _sparse_symbolic_factor(
    ftype::SparseFactorization_t,
    structure::SparseMatrixStructure,
    sym_opts::SparseSymbolicFactorOptions,
)::SparseOpaqueSymbolicFactorization
    return @ccall LIBSPARSE._Z12SparseFactorh21SparseMatrixStructure27SparseSymbolicFactorOptions(
        ftype::Cuint,
        structure::SparseMatrixStructure,
        sym_opts::SparseSymbolicFactorOptions,
    )::SparseOpaqueSymbolicFactorization
end

# Numeric factor on top of an existing symbolic factor. Reusable: the
# symbolic handle is not consumed.
function _sparse_numeric_factor(
    symbolic::SparseOpaqueSymbolicFactorization,
    matrix::SparseMatrix_t,
    num_opts::SparseNumericFactorOptions,
)::SparseOpaqueFactorization_t
    return @ccall LIBSPARSE._Z12SparseFactor33SparseOpaqueSymbolicFactorization19SparseMatrix_Double26SparseNumericFactorOptions(
        symbolic::SparseOpaqueSymbolicFactorization,
        matrix::SparseMatrix_t,
        num_opts::SparseNumericFactorOptions,
    )::SparseOpaqueFactorization_t
end

# In-place solve `A · X = B`, where B is a column-major dense matrix and X
# overwrites B's storage.
function _sparse_solve_matrix!(
    factor::SparseOpaqueFactorization_t,
    B::DenseMatrix_t,
)
    @ccall LIBSPARSE._Z11SparseSolve32SparseOpaqueFactorization_Double18DenseMatrix_Double(
        factor::SparseOpaqueFactorization_t,
        B::DenseMatrix_t,
    )::Cvoid
    return nothing
end

function _sparse_solve_vector!(
    factor::SparseOpaqueFactorization_t,
    b::DenseVector_t,
)
    @ccall LIBSPARSE._Z11SparseSolve32SparseOpaqueFactorization_Double18DenseVector_Double(
        factor::SparseOpaqueFactorization_t,
        b::DenseVector_t,
    )::Cvoid
    return nothing
end

# `Y = A · X`, dense multi-column. `Y` must be allocated to (rowCount, ncols)
# by the caller. libSparse overwrites — does not accumulate.
function _sparse_multiply_matrix!(
    A::SparseMatrix_t,
    X::DenseMatrix_t,
    Y::DenseMatrix_t,
)
    @ccall LIBSPARSE._Z14SparseMultiply19SparseMatrix_Double18DenseMatrix_DoubleS0_(
        A::SparseMatrix_t,
        X::DenseMatrix_t,
        Y::DenseMatrix_t,
    )::Cvoid
    return nothing
end

function _sparse_multiply_vector!(
    A::SparseMatrix_t,
    x::DenseVector_t,
    y::DenseVector_t,
)
    @ccall LIBSPARSE._Z14SparseMultiply19SparseMatrix_Double18DenseVector_DoubleS0_(
        A::SparseMatrix_t,
        x::DenseVector_t,
        y::DenseVector_t,
    )::Cvoid
    return nothing
end

# Frees the libSparse-side numeric / symbolic storage attached to an opaque
# factor. Idempotent: a second call with a `SparseStatusReleased` handle is
# a no-op on libSparse's side.
function _sparse_cleanup_factor!(factor::SparseOpaqueFactorization_t)
    @ccall LIBSPARSE._Z13SparseCleanup32SparseOpaqueFactorization_Double(
        factor::SparseOpaqueFactorization_t,
    )::Cvoid
    return nothing
end

function _sparse_cleanup_symbolic!(symbolic::SparseOpaqueSymbolicFactorization)
    @ccall LIBSPARSE._Z13SparseCleanup33SparseOpaqueSymbolicFactorization(
        symbolic::SparseOpaqueSymbolicFactorization,
    )::Cvoid
    return nothing
end

# Translate libSparse status codes into Julia exceptions. Singular and
# parameter-error are the most common; the rest fall through to a generic
# `error`.
function _libsparse_throw(status::SparseStatus_t, op::AbstractString)
    status == SparseMatrixIsSingular &&
        throw(LinearAlgebra.SingularException(0))
    status == SparseParameterError &&
        throw(ArgumentError("libSparse $(op) failed: parameter error"))
    status == SparseInternalError &&
        error("libSparse $(op) failed: internal error")
    status == SparseStatusFailed &&
        error("libSparse $(op) failed")
    return error("libSparse $(op) failed: status=$(Int(status))")
end
