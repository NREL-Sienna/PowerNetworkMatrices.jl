"""
Line Outage Distribution Factors (LODFs) are a sensitivity measure of how a change in
a line’s flow affects the flows on other lines in the system.
"""
struct LODF{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Float64}
    data::Array{Float64, 2}
    axes::Ax
    lookup::L
end

function _buildlodf(branches, nodes, bus_lookup::Dict{Int, Int}, dist_slack::Array{Float64} = [0.1])
    linecount = length(branches)
    ptdf, a = calculate_PTDF_matrix_DENSE(branches, nodes, bus_lookup, dist_slack)
    @error size(a)
    ptdf_denominator = ptdf*a
    #ptdf_denominator = gemm('N', 'N', ptdf, Matrix{Float64}(a))
    for iline in 1:linecount
        if (1.0 - ptdf_denominator[iline, iline]) < 1.0E-06
            ptdf_denominator[iline, iline] = 0.0
        end
    end
    (Dem, dipiv, dinfo) = getrf!(
        Matrix{Float64}(LinearAlgebra.I, linecount, linecount) -
        Array(LinearAlgebra.Diagonal(ptdf_denominator)),
    )
    lodf = gemm('N', 'N', H, getri!(Dem, dipiv))
    lodf =
        lodf - Array(LinearAlgebra.Diagonal(lodf)) -
        Matrix{Float64}(LinearAlgebra.I, linecount, linecount)
    return lodf
end

"""
Builds the LODF matrix from a group of branches and nodes. The return is a LOLDF array indexed with the branch name.

# Keyword arguments
- `dist_slack::Vector{Float64}`: Vector of weights to be used as distributed slack bus.
    The distributed slack vector has to be the same length as the number of buses
"""
function LODF(branches, nodes, dist_slack::Vector{Float64} = [0.1])

    #Get axis names
    line_ax = [branch.name for branch in branches]
    axes = (line_ax, line_ax)
    look_up = (make_ax_ref(line_ax), make_ax_ref(line_ax))
    bus_ax = [PSY.get_number(bus) for bus in nodes]
    lodf = _buildlodf(branches, nodes, make_ax_ref(bus_ax), dist_slack)
    return LODF(lodf, axes, look_up)
end

"""
Builds the LODF matrix from a system. The return is a LOLDF array indexed with the branch name.

# Keyword arguments
- `dist_slack::Vector{Float64}`: Vector of weights to be used as distributed slack bus.
    The distributed slack vector has to be the same length as the number of buses
"""
function LODF(sys::PSY.System, dist_slack::Vector{Float64} = [0.1])
    branches = sort!(
        collect(PSY.get_components(PSY.ACBranch, sys));
        by = x ->
            (PSY.get_number(PSY.get_arc(x).from), PSY.get_number(PSY.get_arc(x).to)),
    )
    nodes = sort!(collect(PSY.get_components(PSY.Bus, sys)); by = x -> PSY.get_number(x))
    return LODF(branches, nodes, dist_slack)
end