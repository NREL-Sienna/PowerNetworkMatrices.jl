# BA matrix ##################################################################
struct BA_Matrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Float64}
    data::SparseArrays.SparseMatrixCSC{Float64, Int}
    axes::Ax
    lookup::L
    ref_bus_positions::Vector{Int}
end

# create the BA matrix
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

# ABA matrix #################################################################
struct ABA_Matrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Float64}
    data::SparseArrays.SparseMatrixCSC{Float64, Int}
    axes::Ax
    lookup::L
    ref_bus_positions::Vector{Int}
end

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
