# BA matrix ##################################################################
struct BA_Matrix <: PowerNetworkMatrix{Float64}
    data::SparseArrays.SparseMatrixCSC{Float64, Int32}
end

# create the BA matrix
function BA_Matrix(sys::PSY.System)
    branches = get_ac_branches(sys)
    buses = get_buses(sys)
    slack_positions = find_slack_positions(buses)
    bus_lookup = _make_ax_ref(buses)

    data = calculate_BA_matrix(branches, slack_positions, bus_lookup)
    return BA_Matrix(data)
end

# ABA matrix #################################################################
struct ABA_Matrix <: PowerNetworkMatrix{Float64}
    data::SparseArrays.SparseMatrixCSC{Float64, Int32}
end

# create the ABA matrix
function ABA_Matrix(A::SparseArrays.SparseMatrixCSC{Int8, Int32},
    BA::SparseArrays.SparseMatrixCSC{T, Int32}
    where {T <: Union{Float32, Float64}},
    slack_positions::Vector{Int64})
    data = calculate_ABA_matrix(A, BA, slack_positions)
    return ABA_Matrix(data)
end
