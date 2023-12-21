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
- `radial_branches::RadialBranches`:
        Structure containing the radial branches and leaf buses that were removed
        while evaluating the matrix
"""
struct IncidenceMatrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Int8}
    data::SparseArrays.SparseMatrixCSC{Int8, Int}
    axes::Ax
    lookup::L
    ref_bus_positions::Set{Int}
    radial_branches::RadialBranches
end

# functions to get stored data
get_axes(A::IncidenceMatrix) = A.axes
get_lookup(A::IncidenceMatrix) = A.lookup
get_slack_position(A::IncidenceMatrix) = A.ref_bus_positions

"""
Builds the Incidence matrix of the system by evaluating the actual matrix and other relevant
values.

# Arguments
- `sys::PSY.System`: the PowerSystem system to consider
- `radial_branches::RadialBranches`:
        Structure containing the radial branches and leaf buses that were removed
        while evaluating the matrix
"""
function evaluate_A_matrix_values(
    sys::PSY.System,
    radial_branches::RadialBranches=RadialBranches()
)
    branches = get_ac_branches(sys, radial_branches.radial_branches)
    buses = get_buses(sys, radial_branches.bus_reduction_map)
    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in buses]
    data, ref_bus_positions = calculate_A_matrix(branches, buses)
    axes = (line_ax, bus_ax)
    lookup = (make_ax_ref(line_ax), make_ax_ref(bus_ax))
    return data, axes, lookup, ref_bus_positions
end

"""
Builds the Incidence matrix of the system by evaluating the actual matrix and other relevant
values.

# Arguments
- `sys::PSY.System`:
        the PowerSystem system to consider
- `reduce_radial_branches::Bool`:
        if True the matrix will be evaluated discarding
        all the radial branches and leaf buses (optional, default value is false)
"""
function IncidenceMatrix(
    sys::PSY.System;
    reduce_radial_branches::Bool=false
)
    data, axes, lookup, ref_bus_positions = evaluate_A_matrix_values(sys)
    if reduce_radial_branches
        rb = RadialBranches(data, lookup[1], lookup[2], ref_bus_positions)
        data, axes, lookup, ref_bus_positions = evaluate_A_matrix_values(sys, rb)
    else
        rb = RadialBranches()
    end
    return IncidenceMatrix(data, axes, lookup, ref_bus_positions, rb)
end

# this function needs to stay here due to dependency of IncidenceMatrix
function RadialBranches(A::IncidenceMatrix)
    return RadialBranches(A.data, A.lookup[1], A.lookup[2], A.ref_bus_positions)
end

