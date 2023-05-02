"""
Line Outage Distribution Factors (LODFs) are a sensitivity measure of how a change in
a line's flow affects the flows on other lines in the system.

# Arguments
- `data<:AbstractArray{Float64, 2}`:
        the actual Incidence matrix.
- `axes<:NTuple{2, Dict}`:
        Tuple containing two vectors (the first one showing the branches names,
        the second showing the buses numbers).
- `lookup<:NTuple{2, Dict}`:
        Tuple containing two dictionaries, the first mapping the branches
        and buses with their enumerated indexes.
"""
struct LODF{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Float64}
    data::Array{Float64, 2}
    axes::Ax
    lookup::L
end

function _buildlodf(
    a::SparseArrays.SparseMatrixCSC{Int8, Int},
    ptdf::Matrix{Float64},
    linear_solver = "KLU",
)
    if linear_solver == "KLU"
        lodf_t = _calculate_LODF_matrix_KLU(a, ptdf)
    elseif linear_solver == "Dense"
        lodf_t = _calculate_LODF_matrix_DENSE(a, ptdf)
    end

    return lodf_t
end

function _calculate_LODF_matrix_KLU(
    a::SparseArrays.SparseMatrixCSC{Int8, Int},
    ptdf::Matrix{Float64},
)
    linecount = size(ptdf, 2)
    ptdf_denominator_t = a * ptdf
    m_I = Int[]
    m_J = Int[]
    m_V = Float64[]
    for iline in 1:linecount
        if (1.0 - ptdf_denominator_t[iline, iline]) < 1.0E-06
            push!(m_I, iline)
            push!(m_J, iline)
            push!(m_V, 1.0)
        else
            push!(m_I, iline)
            push!(m_J, iline)
            push!(m_V, 1 - ptdf_denominator_t[iline, iline])
        end
    end
    Dem_LU = klu(SparseArrays.sparse(m_I, m_J, m_V))
    lodf_t = Dem_LU \ ptdf_denominator_t
    lodf_t[SparseArrays.diagind(lodf_t)] .= -1

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
Builds the LODF matrix from a group of branches and buses. The return is a LOLDF array indexed with the branch name.

# Arguments
- `branches`:
        vector of the System AC branches
- `buses::Vector{PSY.Bus}`:
        vector of the System buses
- `dist_slack::Vector{Float64}`:
        Vector of weights to be used as distributed slack bus.
        The distributed slack vector has to be the same length as the number of buses.

"""
function LODF(branches, buses::Vector{PSY.Bus}, dist_slack::Vector{Float64} = Float64[])

    # get axis names
    line_ax = [branch.name for branch in branches]
    axes = (line_ax, line_ax)
    look_up = (make_ax_ref(line_ax), make_ax_ref(line_ax))
    bus_ax = [PSY.get_number(bus) for bus in buses]
    lodf = _buildlodf(branches, buses, make_ax_ref(bus_ax), dist_slack)
    return LODF(lodf, axes, look_up)
end

"""
Builds the LODF matrix from a system. The return is a LOLDF array indexed with the branch name.

# Arguments
- `sys::PSY.System`:
        Power Systems system
- `dist_slack::Vector{Float64}`:
        Vector of weights to be used as distributed slack bus.
        The distributed slack vector has to be the same length as the number of buses.
"""
function LODF(
    sys::PSY.System,
    linear_solver::String = "KLU",
    dist_slack::Vector{Float64} = Float64[],
)
    branches = get_ac_branches(sys)
    nodes = get_buses(sys)
    return LODF(branches, nodes, linear_solver, dist_slack)
end

function LODF(
    A::IncidenceMatrix,
    PTDFm::PTDF,
    linear_solver::String = "KLU",
)
    ax_ref = make_ax_ref(A.axes[1])
    return LODF(
        _buildlodf(A.data, PTDFm.data, linear_solver),
        (A.axes[1], A.axes[1]),
        (ax_ref, ax_ref),
    )
end
