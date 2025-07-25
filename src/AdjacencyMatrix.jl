
"""
Nodal incidence matrix (Adjacency) is an N x N matrix describing a power system with N buses. It represents the directed connectivity of the buses in a power system.

The AdjacencyMatrix Struct is indexed using the Bus Numbers, no need for them to be sequential

# Arguments
- `data::SparseArrays.SparseMatrixCSC{Int8, Int}`:
        stores the incidence matrix
- `axes<:NTuple{2, Dict}`:
        Tuple containing two vectors, the first one contains the names of each
        line of the network (each one related to a row of the Matrix in "data"),
        the second one contains the names of each bus of the network (each one
        related to a column of the Matrix in "data")
- `lookup<:NTuple{2, Dict}`:
        Tuple containing 2 Dictionaries mapping the number of rows and columns
        with the names of branches and buses
- `ref_bus_positions::Set{Int}`:
        Vector containing the indexes of the columns of the BA matrix corresponding
        to the reference buses
"""
struct AdjacencyMatrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Int8}
    data::SparseArrays.SparseMatrixCSC{Int8, Int}
    axes::Ax
    lookup::L
    ref_bus_positions::Set{Int}
    network_reduction_data::NetworkReductionData
end

# functions to get stored data
get_axes(A::AdjacencyMatrix) = A.axes
get_lookup(A::AdjacencyMatrix) = A.lookup
get_slack_position(A::AdjacencyMatrix) = A.ref_bus_positions
get_network_reduction_data(A::AdjacencyMatrix) = A.network_reduction_data
get_bus_axis(A::AdjacencyMatrix) = A.axes[1]
get_bus_lookup(A::AdjacencyMatrix) = A.lookup[1]

function get_reduction(
    A::AdjacencyMatrix,
    sys::PSY.System,
    reduction::DegreeTwoReduction,
)
    irreducible_buses = get_irreducible_buses(reduction)
    validate_buses(A, irreducible_buses)
    network_reduction_data = get_network_reduction_data(A)
    direct_branch_map = network_reduction_data.direct_branch_map
    parallel_branch_map = network_reduction_data.parallel_branch_map
    transformer3W_map = network_reduction_data.transformer3W_map

    for c in PSY.get_components(PSY.StaticInjection, sys)
        bus = PSY.get_bus(c)
        if PSY.get_available(bus)
            push!(irreducible_buses, PSY.get_number(bus))
        end
    end
    exempt_bus_positions = Set(get_irreducible_indices(A, irreducible_buses))
    series_branch_map, reverse_series_branch_map, removed_buses, removed_arcs =
        get_degree2_reduction(
            A.data,
            A.lookup[2],
            exempt_bus_positions,
            direct_branch_map,
            parallel_branch_map,
            transformer3W_map,
        )
    return NetworkReductionData(;
        irreducible_buses = Set(irreducible_buses),
        series_branch_map = series_branch_map,
        reverse_series_branch_map = reverse_series_branch_map,
        removed_buses = removed_buses,
        removed_arcs = removed_arcs,
        reductions = ReductionContainer(; degree_two_reduction = reduction),
    )
end

"""
Builds a AdjacencyMatrix from the system. The return is an N x N AdjacencyMatrix Array indexed with the bus numbers.

# Arguments
- `check_connectivity::Bool`:
        Checks connectivity of the network using Depth First Search (DFS) algorithm
"""
function AdjacencyMatrix(sys::PSY.System; check_connectivity::Bool = true, kwargs...)
    ybus = Ybus(
        sys;
        check_connectivity = check_connectivity,
        kwargs...,
    )
    return AdjacencyMatrix(ybus)
end

function AdjacencyMatrix(ybus::Ybus)
    adj_matrix = deepcopy(ybus.adjacency_data)
    for i in 1:size(adj_matrix, 1)
        adj_matrix[i, i] = 0
    end
    SparseArrays.dropzeros!(adj_matrix)
    return AdjacencyMatrix(
        adj_matrix,
        deepcopy(ybus.axes),
        deepcopy(ybus.lookup),
        Set([ybus.lookup[1][x] for x in ybus.ref_bus_numbers]),
        ybus.network_reduction_data,
    )
end

"""
Validates connectivity by checking that the number of subnetworks is 1 (fully connected network).
"""
function validate_connectivity(M::AdjacencyMatrix)
    sub_nets = find_subnetworks(M)
    return length(sub_nets) == 1
end

"""
Evaluates subnetworks by looking for the subsets of buses connected each other,
but not connected with buses of other subsets.
"""
function find_subnetworks(M::AdjacencyMatrix)
    bus_numbers = get_bus_axis(M)
    return find_subnetworks(M.data, bus_numbers)
end
