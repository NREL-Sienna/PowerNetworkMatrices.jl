"""
Incidence matrix: shows connection between buses, defining lines
"""

# define structure for incidence matrix A
struct IncidenceMatrix <: PowerNetworkMatrix{Int}
    data::SparseArrays.SparseMatrixCSC{Int8, Int32}
    axes::Tuple{Vector{Int64}, Vector{Int64}}
    lookup::Dict{Int64, Int64}
    slack_positions::Vector{Int}
end

# functions to get stored data
get_axes(A::IncidenceMatrix) = A.axes
get_lookup(A::IncidenceMatrix) = A.lookup
get_slack_position(A::IncidenceMatrix) = A.slack_positions

# create the incidence matrix A
function IncidenceMatrix(sys::PSY.System)
    branches = get_ac_branches(sys)
    buses = get_buses(sys)
    # ! column related to slack bus has been removed
    data, _ = calculate_A_matrix(branches, buses)
    slack_positions = find_slack_positions(buses)
    lookup = _make_ax_ref(buses)
    axes = (PSY.get_number.(branches), PSY.get_number.(buses))
    return IncidenceMatrix(data, axes, lookup, slack_positions)
end
