function A_mul_B!(
    C::Matrix{Float64},
    A::SparseArrays.SparseMatrixCSC{Float64, Int},
    B::Matrix{Float64},
    ref_bus_positions::Set{Int},
)
    nzv = SparseArrays.nonzeros(A)
    rv = SparseArrays.rowvals(A)
    colptr = SparseArrays.getcolptr(A)
    for k in 1:size(C, 2)
        k in ref_bus_positions && continue
        offset_1 = sum(k .> ref_bus_positions)
        @inbounds for col in 1:size(A, 2)
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

"""
temporary: used for evaluation. Next commits need to chose between this and A_mul_B!
"""

function A_mul_B_2!(
    C::Matrix{Float64},
    A::SparseArrays.SparseMatrixCSC{Float64, Int},
    B::Matrix{Float64},
    ref_bus_positions::Set{Int},
)
    # get sparse matrix values
    nzv = SparseArrays.nonzeros(A)
    rv = SparseArrays.rowvals(A)
    colptr = SparseArrays.getcolptr(A)
    # evaluate the ptdf index positions that are not ref buses
    col_idx = Vector{Int}(undef, size(A, 2) - length(ref_bus_positions))
    sort!(ref_bus_positions)
    offset = 0
    val_1 = 1
    for i in 1:length(col_idx)
        if val_1 <= length(ref_bus_positions)
            if i + offset == ref_bus_positions[1 + offset]
                offset += 1
                val_1 = 1 + offset
            end
        end
        col_idx[i] = i + offset
    end
    # start loop
    for k in 1:length(col_idx)
        for col in 1:length(col_idx)
            B_val = B[col, k]
            for j in colptr[col_idx[col]]:(colptr[col_idx[col] + 1] - 1)
                C[rv[j], col_idx[k]] += nzv[j] * B_val
            end
        end
    end
    return
end
