# BA matrix ##################################################################
struct BA_Matrix <: PowerNetworkMatrix{Float64}
    data::SparseArrays.SparseMatrixCSC{Float64, Int}
    ref_bus_positions::Vector{Int}
end

# create the BA matrix
function BA_Matrix(sys::PSY.System)
    branches = get_ac_branches(sys)
    buses = get_buses(sys)
    ref_bus_positions = find_slack_positions(buses)
    bus_lookup = make_ax_ref(buses)

    data = calculate_BA_matrix(branches, ref_bus_positions, bus_lookup)
    return BA_Matrix(data, ref_bus_positions)
end

# ABA matrix #################################################################
struct ABA_Matrix <: PowerNetworkMatrix{Float64}
    data::SparseArrays.SparseMatrixCSC{Float64, Int}
end

# create the ABA matrix
function ABA_Matrix(A::SparseArrays.SparseMatrixCSC{Int8, Int},
    BA::SparseArrays.SparseMatrixCSC{T, Int},
    slack_positions::Vector{Int}) where {T <: Union{Float32, Float64}}
    data = calculate_ABA_matrix(A, BA, slack_positions)
    return ABA_Matrix(data)
end
