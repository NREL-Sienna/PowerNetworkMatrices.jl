const DC_ABA_Matrix_Factorized = ABA_Matrix{
    Tuple{Vector{Int64}, Vector{Int64}},
    Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
    KLU.KLUFactorization{Float64, Int64},
}
const DC_ABA_Matrix_Unfactorized = ABA_Matrix{
    Tuple{Vector{Int64}, Vector{Int64}},
    Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
    Nothing,
}
# TODO: what about MKLPardiso? switch to <:LinearAlgebra.Factorization instead?
@static if USE_AA
    const DC_vPTDF_Matrix = VirtualPTDF{
        Tuple{Vector{Tuple{Int, Int}}, Vector{Int64}},
        Tuple{Dict{Tuple{Int, Int}, Int64}, Dict{Int64, Int64}},
        <:Union{
            AppleAccelerate.AAFactorization{Float64},
            KLU.KLUFactorization{Float64, Int64},
        },
    }
else
    const DC_vPTDF_Matrix = VirtualPTDF{
        Tuple{Vector{Tuple{Int, Int}}, Vector{Int64}},
        Tuple{Dict{Tuple{Int, Int}, Int64}, Dict{Int64, Int64}},
        KLU.KLUFactorization{Float64, Int64},
    }
end
const DC_PTDF_Matrix = PTDF{
    Tuple{Vector{Int64}, Vector{Tuple{Int, Int}}},
    Tuple{Dict{Int64, Int64}, Dict{Tuple{Int, Int}, Int64}},
    Matrix{Float64},
}
const DC_BA_Matrix = BA_Matrix{
    Tuple{Vector{Int64}, Vector{Tuple{Int, Int}}},
    Tuple{Dict{Int64, Int64}, Dict{Tuple{Int, Int}, Int64}},
}
const AC_Ybus_Matrix = Ybus{
    Tuple{Vector{Int64}, Vector{Int64}},
    Tuple{Dict{Int64, Int64}, Dict{Int64, Int64}},
}
