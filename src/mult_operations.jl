function A_mul_B!(
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
    col_idx = setdiff(1:size(A, 2), ref_bus_positions)
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

function A_div_B!(
    C::Matrix{Float64},
    K::KLU.KLUFactorization{Float64, Int64},
    B::SparseArrays.SparseMatrixCSC{Float64, Int64},
    ref_bus_positions::Set{Int},
)
    # get sparse matrix values
    nzv = SparseArrays.nonzeros(A)
    rv = SparseArrays.rowvals(A)
    colptr = SparseArrays.getcolptr(A)
    # evaluate the ptdf index positions that are not ref buses
    col_idx = setdiff(1:size(A, 2), ref_bus_positions)
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