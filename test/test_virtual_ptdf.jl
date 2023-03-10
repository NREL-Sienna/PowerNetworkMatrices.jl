# include("test_data.jl")

# @testset "Test indexing methods" begin
# end

# @testset "Test actual numbers" begin
# end

# first test on the vPTDF construction

using Revise
using PowerSystems
using PowerSystemCaseBuilder
using PowerNetworkMatrices
using Base


sys = System("ACTIVSg2000.m")

ptdf_complete = PTDF(sys; linear_solver="KLU")

ptdf_virtual = VirtualPTDF(sys)

# # These below might not be needed
Base.getindex(vptdf::VirtualPTDF, idx::CartesianIndex) = vptdf.data[idx]

function Base.getindex(vptdf::VirtualPTDF, row, column)
    # Here is where the method has to implement the logic calculating the column
    # use the indexes to get the proper entry address
    # implement the writing to cache and so on

    # at first check if value is in the cache


    # evaluate the value for the PTDF column

    return
end


