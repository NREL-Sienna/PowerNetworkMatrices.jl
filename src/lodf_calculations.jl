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
    branches,
    buses::Vector{PSY.Bus},
    bus_lookup::Dict{Int, Int},
    dist_slack::Array{Float64},
)
    ref_bus_positions = find_slack_positions(buses)
    ptdf, a = calculate_PTDF_matrix_KLU(
        branches,
        buses,
        bus_lookup,
        ref_bus_positions,
        dist_slack,
    )
    return _buildlodf(a, ptdf)
end

function _buildlodf(a::SparseArrays.SparseMatrixCSC{Int8, Int}, ptdf::Matrix{Float64})
    linecount = size(ptdf, 2)   # ptdf is buses * lines
    ptdf_denominator = transpose(ptdf) * transpose(a)
    for iline in 1:linecount
        if (1.0 - ptdf_denominator[iline, iline]) < 1.0E-06
            ptdf_denominator[iline, iline] = 0.0
        end
    end
    (Dem, dipiv, dinfo) = getrf!(
        Matrix{Float64}(1.0I, linecount, linecount) -
        Array(LinearAlgebra.Diagonal(ptdf_denominator)),
    )
    lodf = gemm('N', 'N', ptdf_denominator, getri!(Dem, dipiv))
    lodf =
        lodf - Array(LinearAlgebra.Diagonal(lodf)) -
        Matrix{Float64}(1.0I, linecount, linecount)
    return lodf
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

    #Get axis names
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
function LODF(sys::PSY.System, dist_slack::Vector{Float64} = Float64[])
    branches = get_ac_branches(sys)
    buses =  get_buses(sys)
    return LODF(branches, buses, dist_slack)
end

function LODF(
    A::IncidenceMatrix,
    PTDFm::PTDF,
)
    ax_ref = make_ax_ref(A.axes[1])
    return LODF(_buildlodf(A.data, PTDFm.data), (A.axes[1], A.axes[1]), (ax_ref, ax_ref))
end
