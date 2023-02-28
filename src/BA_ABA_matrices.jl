# BA matrix ##################################################################
struct BA_Matrix <: PowerNetworkMatrix{Float64}
    data::SparseArrays.SparseMatrixCSC{Float64, Int32}
end

# create the BA matrix
function BA_Matrix(branches, slack_positions::Vector{Int64},
    bus_lookup::Dict{Int64, Int64})
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
    where {T <: Union{Float32, Float64}})
    data = calculate_ABA_matrix(A, BA)
    return ABA_Matrix(data)
end
