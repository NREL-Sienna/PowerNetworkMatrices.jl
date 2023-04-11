function A_mul_B!(
    C::Matrix{Float64},
    A::SparseArrays.SparseMatrixCSC{Float64, Int},
    B::Matrix{Float64},
    ref_bus_positions::Vector{Int},
)
    nzv = SparseArrays.nonzeros(A)
    rv = SparseArrays.rowvals(A)
    colptr = SparseArrays.getcolptr(A)
    for k in 1:size(C, 2)
        k in ref_bus_positions && continue
        offset_1 = sum(k .> ref_bus_positions)
        for col in 1:size(A, 2)
            col in ref_bus_positions && continue
            offset_2 = sum(col .> ref_bus_positions)
            B_val = B[col - offset_2, k - offset_1]
            for j in colptr[col]:(colptr[col + 1] - 1)
                C[rv[j], k] += nzv[j] * B_val
            end
        end
    end
    return
end
