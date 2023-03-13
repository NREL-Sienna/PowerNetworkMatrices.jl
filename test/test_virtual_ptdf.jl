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
using KLU
using SparseArrays
import LinearAlgebra: ldiv!, mul!

sys = System("ACTIVSg2000.m")

# test if adding data works

# TODO: "setindex" must check if the given values are already contained, and to check that dimensions and name are correct (if ∈ axes)
# TODO: put everything in a test file and in the functions directory

# These below might not be needed
function _buildptdf_columns(
    K::KLU.KLUFactorization{Float64, Int32},
    BA_row::SparseArrays.SparseVector{Float64, Int32})
    # inizialize matrices for evaluation
    col = zeros(length(BA_row))
    ldiv!(col, K, BA_row)

    return col
end

# functions to get indeces
function _get_line_index(vptdf::VirtualPTDF, row::String)
    row_ = findall(x->x==row, vptdf.axes[1])
    if length(row_) > 1
        error("multiple lines with the same name $row in vptdf.axes[1]")
    else
        row = row_[1]
    end
    return row
end

function _get_value(vptdf::VirtualPTDF, row::Union{Int, String}, column::Union{Int, Colon})
    # get stored value if present
    if !isempty(vptdf.data[1][row]) && !isempty(vptdf.data[2][row])
        vptdf.data[3][row] =+ 1
        if column isa Colon
            return vptdf.data[2][row]
        else
            return vptdf.data[2][row][column]
        end
    else
        return []
    end
end

function Base.getindex(vptdf::VirtualPTDF, row::Union{Int, String}, column::Union{Int, Colon})
    # at first get the index value if needed
    if row isa String
        row = _get_line_index(vptdf, row)
    end
    # at first check if value is in the cache
    value_ = _get_value(vptdf, row, column)
    if !isempty(value_)
        return value_
    else
        # evaluate the value for the PTDF column
        PTDF_row = _buildptdf_columns(vptdf.K, vptdf.BA[row, :])
        # add slack bus value (zero)
        PTDF_row = vcat(PTDF_row[1:vptdf.slack_positions[1]-1], [0], PTDF_row[vptdf.slack_positions[1]:end])
        # add the valute in cache
        vptdf.data[1][row] = vptdf.axes[1][row]
        vptdf.data[2][row] = PTDF_row
        vptdf.data[3][row] =+ 1
        # return value
        if column isa Colon
            return PTDF_row
        else
            return PTDF_row[column]
        end
    end
end

# functions to set indeces

## tests #####################################################################

ptdf_complete = PTDF(sys; linear_solver="KLU")
ptdf_virtual = VirtualPTDF(sys)

count_ = 0
for i in axes(ptdf_complete.data, 1)
    for j in axes(ptdf_complete.data, 2)
        # check values
        if !isapprox(ptdf_complete.data[i, j], ptdf_virtual[i, j], atol=1e-10)
            @show count_ =+ 1
            # @show ptdf_complete.data[i, j] - ptdf_virtual[i, j]
        end
    end
end

ptdf_virtual = VirtualPTDF(sys)
count1_ = 0
for i in axes(ptdf_complete.data, 1)
    # check values
    if !isapprox(ptdf_complete.data[i, :], ptdf_virtual[i, :], atol=1e-10)
        @show count1_ =+ 1
        # @show ptdf_complete.data[i, j] - ptdf_virtual[i, j]
    end
end

ptdf_virtual = VirtualPTDF(sys)
count2_ = 0
for (i, name) in enumerate(ptdf_complete.axes[1])
    if !isapprox(ptdf_complete.data[i, :], ptdf_virtual[name, :], atol=1e-10)
        @show count2_ =+ 1
        # @show ptdf_complete.data[i, j] - ptdf_virtual[i, j]
    end
end

"""
no differences
"""

# final check on entire matrix
new_ptdf = reduce(vcat,transpose.(ptdf_virtual.data[2]))
old_ptds = ptdf_complete.data
isapprox(new_ptdf, old_ptds, atol=1e-10)