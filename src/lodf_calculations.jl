"""
The Line Outage Distribution Factor (LODF) matrix gathers a sensitivity coefficients
of how a change in a line's flow affects the flows on other lines in the system.

# Arguments
- `data<:AbstractArray{Float64, 2}`:
        the transposed LODF matrix.
- `axes<:NTuple{2, Dict}`:
        Tuple containing two identical vectors containing the names of the 
        branches related to each row/column.
- `lookup<:NTuple{2, Dict}`:
        Tuple containing two identical dictionaries mapping the branches
         their enumerated indexes (row and column numbers).
- `tol::Base.RefValue{Float64}`:
        tolerance used for sparsifying the matrix (dropping element whose
        absolute value is below this threshold).
"""
struct LODF{Ax, L <: NTuple{2, Dict}, M <: AbstractArray{Float64, 2}} <:
       PowerNetworkMatrix{Float64}
    data::M
    axes::Ax
    lookup::L
    tol::Base.RefValue{Float64}
end

"""
Sets to zero those elements of the LODF matrix whose absolute values are below
the threshold specified by the field "tol".

# Arguments
- `mat::LODF`:
        LODF structure
- `tol::Float64`:
        tolerance
"""
function drop_small_entries!(mat::LODF, tol::Float64)
    if tol < mat.tol[]
        @info "Specified tolerance is smaller than the current tolerance."
    end
    make_entries_zero!(mat.data, tol)
    mat.tol[] = tol
    return
end

function _buildlodf(
    a::SparseArrays.SparseMatrixCSC{Int8, Int},
    ptdf::Matrix{Float64},
    linear_solver::String = "KLU",
)
    if linear_solver == "KLU"
        lodf_t = _calculate_LODF_matrix_KLU(a, ptdf)
    elseif linear_solver == "Dense"
        lodf_t = _calculate_LODF_matrix_DENSE(a, ptdf)
    else
        error("MKLPardiso still to be implemented")
    end
    return lodf_t
end

function _buildlodf(
    a::SparseArrays.SparseMatrixCSC{Int8, Int64},
    k::KLU.KLUFactorization{Float64, Int64},
    ba::SparseArrays.SparseMatrixCSC{Float64, Int64},
    ref_bus_positions::Set{Int},
    linear_solver::String,
)
    if linear_solver == "KLU"
        lodf_t = _calculate_LODF_matrix_KLU(a, k, ba, ref_bus_positions)
    else
        error("Other methods still to be implemented.")
    end
    return lodf_t
end

function _calculate_LODF_matrix_KLU(
    a::SparseArrays.SparseMatrixCSC{Int8, Int64},
    k::KLU.KLUFactorization{Float64, Int64},
    ba::SparseArrays.SparseMatrixCSC{Float64, Int64},
    ref_bus_positions::Set{Int64},
)
    linecount = size(ba, 2)
    # get inverse of aba
    first_ = zeros(size(a, 2), size(a, 1))
    valid_ix = setdiff(1:size(a, 2), ref_bus_positions)
    copyto!(first_, transpose(a))
    first_[valid_ix, :] = KLU.solve!(k, first_[valid_ix, :])
    first_[collect(ref_bus_positions), :] .= 0.0
    ptdf_denominator = first_' * ba

    m_I = Int[]
    m_V = Float64[]
    for iline in 1:linecount
        if (1.0 - ptdf_denominator[iline, iline]) < 1.0E-06
            push!(m_I, iline)
            push!(m_V, 1.0)
        else
            push!(m_I, iline)
            push!(m_V, 1 - ptdf_denominator[iline, iline])
        end
    end
    Dem_LU = klu(SparseArrays.sparse(m_I, m_I, m_V))
    lodf_t = Dem_LU \ ptdf_denominator
    lodf_t[SparseArrays.diagind(lodf_t)] .= -1
    return lodf_t
end

function _calculate_LODF_matrix_KLU(
    a::SparseArrays.SparseMatrixCSC{Int8, Int},
    ptdf::Matrix{Float64},
)
    linecount = size(ptdf, 2)
    ptdf_denominator_t = a * ptdf
    m_I = Int[]
    m_V = Float64[]
    for iline in 1:linecount
        if (1.0 - ptdf_denominator_t[iline, iline]) < 1.0E-06
            push!(m_I, iline)
            push!(m_V, 1.0)
        else
            push!(m_I, iline)
            push!(m_V, 1 - ptdf_denominator_t[iline, iline])
        end
    end
    Dem_LU = klu(SparseArrays.sparse(m_I, m_I, m_V))
    lodf_t = Dem_LU \ ptdf_denominator_t
    lodf_t[SparseArrays.diagind(lodf_t)] .= -1

    return lodf_t
end

# ! temp for evaluation
function _calculate_LODF_matrix_KLU2(
    a::SparseArrays.SparseMatrixCSC{Int8, Int},
    ptdf::Matrix{Float64},
)
    ptdf_denominator_t = a * ptdf
    linecount = size(ptdf, 2)
    lodf_t = zeros(linecount, linecount)
    for i in 1:linecount
        for j in 1:linecount
            if i == j
                lodf_t[i, i] = -1.0
            else
                if (1.0 - ptdf_denominator_t[i, i]) < 1.0E-06
                    lodf_t[i, j] = ptdf_denominator_t[i, j]
                else
                    lodf_t[i, j] = ptdf_denominator_t[i, j] / (1 - ptdf_denominator_t[i, i])
                end
            end
        end
    end
    return lodf_t
end

function _calculate_LODF_matrix_DENSE(
    a::SparseArrays.SparseMatrixCSC{Int8, Int},
    ptdf::Matrix{Float64},
)
    linecount = size(ptdf, 2)
    ptdf_denominator_t = a * ptdf
    for iline in 1:linecount
        if (1.0 - ptdf_denominator_t[iline, iline]) < 1.0E-06
            ptdf_denominator_t[iline, iline] = 0.0
        end
    end
    (Dem, dipiv, dinfo) = getrf!(
        Matrix{Float64}(1.0I, linecount, linecount) -
        Array(LinearAlgebra.Diagonal(ptdf_denominator_t)),
    )
    lodf_t = gemm('N', 'N', getri!(Dem, dipiv), ptdf_denominator_t)
    lodf_t =
        lodf_t - Array(LinearAlgebra.Diagonal(lodf_t)) -
        Matrix{Float64}(1.0I, linecount, linecount)

    return lodf_t
end

"""
Builds the LODF matrix given the data of branches and buses of the system.

# Arguments
- `branches`:
        vector of the System AC branches
- `buses::Vector{PSY.Bus}`:
        vector of the System buses
- `tol::Float64`:
        Tolerance to eliminate entries in the LODF matrix (default eps())
"""
function LODF(
    branches,
    buses::Vector{PSY.Bus};
    linear_solver::String = "KLU",
    tol::Float64 = eps(),
)

    # get axis names
    line_ax = [branch.name for branch in branches]
    axes = (line_ax, line_ax)
    look_up = (make_ax_ref(line_ax), make_ax_ref(line_ax))
    bus_ax = [PSY.get_number(bus) for bus in buses]
    bus_lookup = make_ax_ref(bus_ax)
    # get network matrices
    ptdf_t, a = calculate_PTDF_matrix_KLU(branches, buses, bus_lookup, Float64[])
    if tol > eps()
        lodf_t = _buildlodf(a, ptdf_t, linear_solver)
        return LODF(sparsify(lodf_t, tol), axes, look_up, Ref(tol))
    else
        return LODF(
            _buildlodf(a, ptdf_t, linear_solver),
            axes,
            look_up,
            Ref(tol),
        )
    end
    return
end

"""
Builds the LODF matrix from a system.

# Arguments
- `sys::PSY.System`:
        Power Systems system
"""
function LODF(
    sys::PSY.System;
    kwargs...,
)
    branches = get_ac_branches(sys)
    buses = get_buses(sys)
    return LODF(branches, buses; kwargs...)
end

"""
Builds the LODF matrix given the Incidence Matrix and the PTDF matrix of the system.

NOTE: tol is referred to the LODF sparsification, not the PTDF one. PTDF matrix 
must be considered as NON sparsified ("tol" argument not specified when calling 
the PTDF method).

# Arguments
- `A::IncidenceMatrix`:
        Structure containing the Incidence matrix of the system.
- `PTDFm::PTDF`:
        Strucutre containing the transposed PTDF matrix of the system.
- `linear_solver::String`:
        Linear solver to be used. Options are "Dense" and "KLU".
- `tol::Float64`:
        Tolerance to eliminate entries in the LODF matrix (default eps()).
"""
function LODF(
    A::IncidenceMatrix,
    PTDFm::PTDF;
    linear_solver::String = "KLU",
    tol::Float64 = eps(),
)
    validate_linear_solver(linear_solver)

    if PTDFm.tol.x > 1e-15
        err_msg = string(
            "The argument `tol` in the PTDF matrix was set to a value dirrent than the default one.\n",
            "The PTDF matrix used as input of the LODF matrix must have the default `tol` value.\n",
        )
        error(err_msg)
    end

    ax_ref = make_ax_ref(A.axes[1])
    if tol > eps()
        lodf_t = _buildlodf(A.data, PTDFm.data, linear_solver)
        return LODF(
            sparsify(lodf_t, tol),
            (A.axes[1], A.axes[1]),
            (ax_ref, ax_ref),
            Ref(tol),
        )
    end
    return LODF(
        _buildlodf(A.data, PTDFm.data, linear_solver),
        (A.axes[1], A.axes[1]),
        (ax_ref, ax_ref),
        Ref(tol),
    )
end

"""
Builds the LODF matrix given the Incidence Matrix and the PTDF matrix of the system.

NOTE: this method does not support distributed slack bus.

# Arguments
- `A::IncidenceMatrix`:
        Structure containing the Incidence matrix of the system.
- `ABA::ABA_Matrix`:
        Structure containing the ABA matrix of the system.
- `BA::BA_Matrix`:
        Structure containing the transposed BA matrix of the system.
- `linear_solver::String`:
        Linear solver to be used. Options are "Dense" and "KLU".
- `tol::Float64`:
        Tolerance to eliminate entries in the LODF matrix (default eps()).
"""
function LODF(
    A::IncidenceMatrix,
    ABA::ABA_Matrix,
    BA::BA_Matrix;
    linear_solver::String = "KLU",
    tol::Float64 = eps(),
)
    validate_linear_solver(linear_solver)
    ax_ref = make_ax_ref(A.axes[1])
    if tol > eps()
        lodf_t = _buildlodf(A.data, ABA.K, BA.data,
            A.ref_bus_positions, linear_solver)
        return LODF(
            sparsify(lodf_t, tol),
            (A.axes[1], A.axes[1]),
            (ax_ref, ax_ref),
            Ref(tol),
        )
    end
    return LODF(
        _buildlodf(A.data, ABA.K, BA.data, A.ref_bus_positions, linear_solver),
        (A.axes[1], A.axes[1]),
        (ax_ref, ax_ref),
        Ref(tol),
    )
end

############################################################
# auxiliary functions for getting data from LODF structure #
############################################################

# NOTE: the LODF matrix is saved as transposed!

function Base.getindex(A::LODF, selected_line, outage_line)
    i, j = to_index(A, outage_line, selected_line)
    return A.data[i, j]
end

function Base.getindex(
    A::LODF,
    selected_line_number::Union{Int, Colon},
    outage_line_number::Union{Int, Colon},
)
    return A.data[outage_line_number, selected_line_number]
end

function get_lodf_data(lodf::LODF)
    return transpose(lodf.data)
end

function get_branch_ax(lodf::LODF)
    return lodf.axes[1]
end

function get_tol(lodf::LODF)
    return lodf.tol
end