"""
Power Transfer Distribution Factors (PTDF) indicate the incremental change in real power that occurs on transmission lines due to real power injections changes at the buses.

The PTDF struct is indexed using the Bus numbers and branch names
"""
struct VirtualPTDF{Ax, L <: NTuple{2, Dict}, T <: Real} <: PowerNetworkMatrix{Real}
    K::KLU.KLUFactorization{Float64, Int32}
    BA::SparseArrays.SparseMatrixCSC{T, Int32}
    slack_positions::Vector{Int64}
    dist_slack::Vector{Float64}
    axes::Ax
    lookup::L
    data::Vector{Vector}

end

function _build_virtualptdf(
    branches,
    nodes::Vector{PSY.Bus})
    slack_positions = find_slack_positions(nodes)
    A, _ = calculate_A_matrix(branches, nodes)
    BA = calculate_BA_matrix(branches, slack_positions, _make_ax_ref(nodes))
    ABA = calculate_ABA_matrix(A, BA, slack_positions)
    K = klu(ABA)

    return BA, K , A, slack_positions

end

"""
Builds the PTDF matrix from a group of branches and nodes. The return is a PTDF array indexed with the bus numbers.

# Keyword arguments
- `dist_slack::Vector{Float64}`: Vector of weights to be used as distributed slack bus.
    The distributed slack vector has to be the same length as the number of buses
"""
function VirtualPTDF(
    branches,
    nodes::Vector{PSY.Bus},
    dist_slack::Vector{Float64} = [0.1])

    #Get axis names
    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in nodes]
    axes = (line_ax, bus_ax)
    look_up = (_make_ax_ref(line_ax), _make_ax_ref(bus_ax))
    BA, K, _, slack_positions = _build_virtualptdf(branches, nodes)
    empty_cache = [["" for i in 1:size(BA,1)], [Float64[] for i in 1:size(BA,1)], zeros(Int32, size(BA,1))]
    return VirtualPTDF(K, BA, slack_positions, dist_slack, axes, look_up, empty_cache)

end

"""
Builds the PTDF matrix from a system. The return is a PTDF array indexed with the bus numbers.

# Keyword arguments
- `dist_slack::Vector{Float64}`: Vector of weights to be used as distributed slack bus.
    The distributed slack vector has to be the same length as the number of buses
"""
function VirtualPTDF(
    sys::PSY.System,
    dist_slack::Vector{Float64} = [0.1])
    branches = get_ac_branches(sys)
    nodes = get_buses(sys)

    return VirtualPTDF(branches, nodes, dist_slack)
end

# overload Base functions
function Base.isempty(vptdf::VirtualPTDF)
    if isempty(vptdf.K.L)
        @warn "Missing L factorization matrix. Use klu(ABA) or klu(A'*BA)."
        return false
    end
    if isempty(vptdf.K.U)
        @warn "Missing U factorization matrix. Use klu(ABA) or klu(A'*BA)."
        return false
    end
    if isempty(vptdf.BA)
        @warn "Missing AB matrix. Use calculate_BA_matrix()."
        return false
    end
    if isempty(vptdf.dist_slack)
        @warn "dist_slack vector is missing."
        return false
    end
    if isempty(vptdf.axes)
        @warn "axes is missing."
        return false
    end
    if isempty(vptdf.lookup)
        @warn "look_up is missing."
        return false
    end
    if isempty(vptdf.data)
        @info "No PTDF rows are stored."
        return false
    end

    return true
end

# Size related overload will work on the BA matrix
Base.size(vptdf::VirtualPTDF) = size(vptdf.BA)
Base.eachindex(vptdf::VirtualPTDF) = CartesianIndices(size(vptdf.BA))

# get data from indeces or line name

Base.getindex(vptdf::VirtualPTDF, idx::CartesianIndex) = vptdf.data[idx]
# Base.getindex(vptdf::VirtualPTDF, rowname::String, bus::Int) = vptdf.data[rowname, bus]

function Base.getindex(vptdf::VirtualPTDF, row, column)
    # Here is where the method has to implement the logic calculating the column
    # use the indexes to get the proper entry address
    # implement the writing to cache and so on

    # at first check if value is in the cache


    return
end



Base.setindex!(A::VirtualPTDF, v, idx...) = A.data[to_index(A, idx...)...] = v
Base.setindex!(A::VirtualPTDF, v, idx::CartesianIndex) = A.data[idx] = v

""" PTDF data is stored in the the cache """
get_data(mat::PowerNetworkMatrix) = mat.data