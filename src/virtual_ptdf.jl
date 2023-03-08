"""
Power Transfer Distribution Factors (PTDF) indicate the incremental change in real power that occurs on transmission lines due to real power injections changes at the buses.

The PTDF struct is indexed using the Bus numbers and branch names
"""
struct VirtualPTDF{Ax, L <: NTuple{2, Dict}, T <: Real} <: PowerNetworkMatrix{Real}
    L::SparseArrays.SparseMatrixCSC{T, Int32}
    U::SparseArrays.SparseMatrixCSC{T, Int32}
    BA::SparseArrays.SparseMatrixCSC{T, Int32}
    axes::Ax
    lookup::L
    cache::Dict
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
    L, U = get_lu_decomposition()
    BA = get_ba_matrix()
    return VirtualPTDF(L, U, BA, axes, look_up, Dict())
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

function Base.isempty(A::VirtualPTDF)
    !isempty(A.L) && return false
    !isempty(A.U) && return false
    !isempty(A.BA) && return false
    !isempty(A.cache) && return false
    return true
end

# Size related overload will work on the BA matrix
Base.size(A::VirtualPTDF) = size(A.BA)
Base.eachindex(A::VirtualPTDF) = CartesianIndices(size(A.BA))

function Base.getindex(A::VirtualPTDF, row, column)
    # Here is where the method has to implement the logic calculating the column
    # use the indexes to get the proper entry address
    # implement the writing to cache and so on
    return
end

# These below might not be needed
Base.getindex(A::VirtualPTDF, idx::CartesianIndex) = A.data[idx]
Base.setindex!(A::VirtualPTDF, v, idx...) = A.data[to_index(A, idx...)...] = v
Base.setindex!(A::VirtualPTDF, v, idx::CartesianIndex) = A.data[idx] = v
