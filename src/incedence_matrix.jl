"""
Incidence matrix: shows connection between buses, defining lines
"""

# define structure for incidence matrix A
struct IncidenceMatrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Int8}
    data::SparseArrays.SparseMatrixCSC{Int8, Int32}
    axes::Ax
    lookup::L
    slack_positions::Vector{Int64}
end

# functions to get stored data
get_axes(A::IncidenceMatrix) = A.axes
get_lookup(A::IncidenceMatrix) = A.lookup
get_slack_position(A::IncidenceMatrix) = A.slack_positions

# create the incidence matrix A
function IncidenceMatrix(sys::PSY.System)
    branches = get_ac_branches(sys)
    buses = get_buses(sys)
    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in buses]
    data, slack_positions = calculate_A_matrix(branches, buses)
    axes = (line_ax, bus_ax)
    lookup = (make_ax_ref(line_ax), make_ax_ref(bus_ax))
    return IncidenceMatrix(data, axes, lookup, slack_positions)
end
