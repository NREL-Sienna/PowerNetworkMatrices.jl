function A_mul_B!(
    C::Matrix{Float64},
    A::SparseArrays.SparseMatrixCSC{Float64, Int},
    B::Matrix{Float64},
    ref_bus_positions::Vector{Int},
)
    for n in 1:size(C, 2)
        n in ref_bus_positions && continue
        offset_1 = sum(n .> ref_bus_positions)
        @inbounds for m in 1:size(A, 2)
            offset_2 = sum(m .> ref_bus_positions)
            B_val = B[m - offset_2, n - offset_1]
            for k in A.colptr[m]:(A.colptr[m + 1] - 1)
                C[A.rowval[k], n] += A.nzval[k] * B_val
            end
        end
    end
    return
end
