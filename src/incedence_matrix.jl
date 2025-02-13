"""
Incidence matrix: shows connection between buses, defining lines

# Arguments
- `data::SparseArrays.SparseMatrixCSC{Int8, Int}`:
        the actual Incidence matrix.
- `axes<:NTuple{2, Dict}`:
        Tuple containing two vectors (the first one showing the branches names,
        the second showing the buses numbers).
- `lookup<:NTuple{2, Dict}`:
        Tuple containing two dictionaries, the first mapping the branches
        and buses with their enumerated indexes.
- `ref_bus_positions::Set{Int}`:
        Vector containing the indices of the reference slack buses.
- `network_reduction::NetworkReduction`:
        Structure containing the details of the network reduction applied when computing the matrix
"""
struct IncidenceMatrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Int8}
    data::SparseArrays.SparseMatrixCSC{Int8, Int}
    axes::Ax
    lookup::L
    ref_bus_positions::Set{Int}
    network_reduction::NetworkReduction
end

# functions to get stored data
get_axes(A::IncidenceMatrix) = A.axes
get_lookup(A::IncidenceMatrix) = A.lookup
get_slack_position(A::IncidenceMatrix) = A.ref_bus_positions

"""
Builds the Incidence matrix of the system by evaluating the actual matrix and other relevant
values.

# Arguments
- `sys::PSY.System`:
        the PowerSystem system to consider
- `network_reduction::NetworkReduction`:
        Structure containing the details of the network reduction applied when computing the matrix
"""
function IncidenceMatrix(
    sys::PSY.System;
    network_reduction::NetworkReduction = NetworkReduction(),
)
    data, axes, lookup, ref_bus_positions = evaluate_A_matrix_values(sys, network_reduction)
    return IncidenceMatrix(data, axes, lookup, ref_bus_positions, network_reduction)
end

"""
Builds the Incidence matrix of the system by evaluating the actual matrix and other relevant
values.

# Arguments
- `sys::PSY.System`: the PowerSystem system to consider
- `network_reduction::NetworkReduction`:
        Structure containing the details of the network reduction applied when computing the matrix
"""
function evaluate_A_matrix_values(
    sys::PSY.System,
    nr::NetworkReduction,
)
    branches = get_ac_branches(sys, nr.removed_branches)
    if !isempty(nr.added_branches)
        branches = vcat(branches, nr.added_branches)
    end
    buses = get_buses(sys, nr.bus_reduction_map)
    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in buses]
    data, ref_bus_positions = calculate_A_matrix(branches, buses)
    axes = (line_ax, bus_ax)
    lookup = (make_ax_ref(line_ax), make_ax_ref(bus_ax))
    return data, axes, lookup, ref_bus_positions
end
