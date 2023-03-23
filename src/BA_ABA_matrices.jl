"""
Structure containing the BA matrix andother relevant data.
Fields:
- data: the BA matrix cominfrom the product between the Incidence Matrix A and
        the Matrix of Susceptance B
- axes: Tuple containing two vectors, the first one contains the names of each 
        line of the network (each one related to a row of the Matrix in "data"),
        the second one contains the names of each bus of the network (each one
        related to a column of the Matrix in "data")
- lookup: 
        Tuple containing 2 Dictionaries mapping the number of rows and columns 
        with the names of branches and buses
- ref_bus_positions:
        Vector containing the indexes of the columns of the BA matrix corresponding
        to the refence buses
"""
struct BA_Matrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Float64}
    data::SparseArrays.SparseMatrixCSC{Float64, Int}
    axes::Ax
    lookup::L
    ref_bus_positions::Vector{Int}
end

"""
Build the BA matrix from a given System
"""
function BA_Matrix(sys::PSY.System)
    branches = get_ac_branches(sys)
    buses = get_buses(sys)
    ref_bus_positions = find_slack_positions(buses)
    bus_lookup = make_ax_ref(buses)
    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in setdiff(buses, ref_bus_positions)]
    axes = (line_ax, bus_ax)
    lookup = (make_ax_ref(line_ax), make_ax_ref(bus_ax))
    data = calculate_BA_matrix(branches, ref_bus_positions, bus_lookup)
    return BA_Matrix(data, axes, lookup, ref_bus_positions)
end

"""
Structure containing the ABA matrix and other relevant data.
Fields:
- data: the ABA matrix cominfrom the product between the Incidence Matrix A and
        the Matrix BA
- axes: Tuple containing two vectors, the first one contains the names of each 
        line of the network (each one related to a row of the Matrix in "data"),
        the second one contains the names of each bus of the network (each one
        related to a column of the Matrix in "data")
- lookup: 
        Tuple containing 2 Dictionaries mapping the number of rows and columns 
        with the names of branches and buses
- ref_bus_positions:
        Vector containing the indexes of the columns of the ABA matrix corresponding
        to the refence buses
"""
struct ABA_Matrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Float64}
    data::SparseArrays.SparseMatrixCSC{Float64, Int}
    axes::Ax
    lookup::L
    ref_bus_positions::Vector{Int}
end

"""
Builds the ABA matrix from a System
"""
function ABA_Matrix(sys::PSY.System)
    branches = get_ac_branches(sys)
    buses = get_buses(sys)
    bus_lookup = make_ax_ref(buses)

    A, ref_bus_positions = calculate_A_matrix(branches, buses)
    BA_full = calculate_BA_matrix_full(branches, bus_lookup)
    data = calculate_ABA_matrix_full(A, BA_full)

    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in buses]
    axes = (line_ax, setdiff(bus_ax, ref_bus_positions))
    lookup = (make_ax_ref(line_ax), make_ax_ref(bus_ax))

    return ABA_Matrix(
        data,
        axes,
        lookup,
        ref_bus_positions,
    )
end
