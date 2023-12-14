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
        vector containing the indices of the reference slack buses.
"""
struct IncidenceMatrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Int8}
    data::SparseArrays.SparseMatrixCSC{Int8, Int}
    axes::Ax
    lookup::L
    ref_bus_positions::Set{Int}
    reduce_radial_branches::RadialBranches
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
- `reduce_radial_branches::RadialBranches`: user defined structure containing the radial
        branches and leaf buses of the system (optional, deafault is empty structure)
"""
function IncidenceMatrix(
    sys::PSY.System,
    radial_branches::Set{String}=Set{String}(),
    bus_reduction_map::Dict{Int, Set{Int}}=Dict{Int, Set{Int}}(),
)
    branches = get_ac_branches(sys, radial_branches)
    buses = get_buses(sys, bus_reduction_map)
    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in buses]
    data, ref_bus_positions = calculate_A_matrix(branches, buses)
    axes = (line_ax, bus_ax)
    lookup = (make_ax_ref(line_ax), make_ax_ref(bus_ax))
    return IncidenceMatrix(data, axes, lookup, ref_bus_positions, RadialBranches(bus_reduction_map, radial_branches))
end

"""
Builds the Incidence matrix of the system by evaluating the actual matrix and other relevant
values.

# Arguments
- `sys::PSY.System`: the PowerSystem system to consider
- `reduce_radial_branches::Bool`: if True the matrix will be evaluated discarding
        all the radial branches and leaf buses (optional, default value is false)
"""
function IncidenceMatrix(
    sys::PSY.System;
    reduce_radial_branches::Bool=false
)
    if reduce_radial_branches
        rb = RadialBranches(IncidenceMatrix(sys))
    else
        rb = RadialBranches()
    end
    return IncidenceMatrix(sys, rb.radial_branches, rb.bus_reduction_map)
end

# this function needs to stay here due to dependency of IncidenceMatrix
function RadialBranches(A::IncidenceMatrix)
    return RadialBranches(A.data, A.lookup[1], A.lookup[2], A.ref_bus_positions)
end

