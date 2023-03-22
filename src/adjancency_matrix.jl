
"""
Nodal incidence matrix (Adjacency) is an N x N matrix describing a power system with N buses. It represents the directed connectivity of the buses in a power system.

The AdjacencyMatrix Struct is indexed using the Bus Numbers, no need for them to be sequential
"""
struct AdjacencyMatrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Int8}
    data::SparseArrays.SparseMatrixCSC{Int8, Int}
    axes::Ax
    lookup::L
    ref_bus_positions::Vector{Int}
end

# functions to get stored data
get_axes(A::AdjacencyMatrix) = A.axes
get_lookup(A::AdjacencyMatrix) = A.lookup
get_slack_position(A::AdjacencyMatrix) = A.ref_bus_positions

"""
Builds a AdjacencyMatrix from the system. The return is an N x N AdjacencyMatrix Array indexed with the bus numbers.

# Keyword arguments
- `check_connectivity::Bool`: Checks connectivity of the network using Goderya's algorithm
"""
function AdjacencyMatrix(sys::PSY.System; check_connectivity::Bool = true, kwargs...)
    nodes = sort!(
        collect(
            PSY.get_components(x -> PSY.get_bustype(x) != BusTypes.ISOLATED, PSY.Bus, sys),
        );
        by = x -> PSY.get_number(x),
    )
    branches = get_ac_branches(sys)
    return AdjacencyMatrix(
        branches,
        nodes;
        check_connectivity = check_connectivity,
        kwargs...,
    )
end

"""
Builds a AdjacencyMatrix from a collection of buses and branches. The return is an N x N AdjacencyMatrix Array indexed with the bus numbers.

# Keyword arguments
- `check_connectivity::Bool`: Checks connectivity of the network using Goderya's algorithm
"""
function AdjacencyMatrix(
    branches,
    nodes::Vector{PSY.Bus};
    check_connectivity::Bool = true,
    kwargs...,
)
    M, bus_lookup = calculate_adjacency(branches, nodes)
    bus_ax = PSY.get_number.(nodes)
    axes = (bus_ax, bus_ax)
    look_up = (bus_lookup, bus_lookup)

    if check_connectivity
        sub_nets = find_subnetworks(M, bus_ax)
        length(sub_nets) > 1 && throw(IS.DataFormatError("Network not connected"))
    end

    return AdjacencyMatrix(M, axes, look_up, find_slack_positions(nodes))
end

function validate_connectivity(M::AdjacencyMatrix)
    sub_nets = find_subnetworks(M)
    return length(sub_nets) == 1
end

function find_subnetworks(M::AdjacencyMatrix)
    bus_numbers = M.axes[2]
    return find_subnetworks(M.data, bus_numbers)
end
