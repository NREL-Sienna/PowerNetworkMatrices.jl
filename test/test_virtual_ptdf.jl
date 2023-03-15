@testset "Virtual PTDF matrices" begin
    sys = PSB.build_system(PSB.PSYTestSystems, "tamu_ACTIVSg2000_sys")
    ptdf_complete = PTDF(sys; linear_solver = "KLU")
    ptdf_virtual = VirtualPTDF(sys)

    for i in axes(ptdf_complete, 1)
        comp = ptdf_complete[i, :]
        virtual = ptdf_virtual[i, :]
        for j in axes(ptdf_complete, 2)
            # check values using PTDFs axes
            @test isapprox(ptdf_complete[i, j], ptdf_virtual[i, j]; atol = 1e-10)
        end
    end
    # Check the cache is populated
    @test length(ptdf_virtual.cache) == length(ptdf_virtual.axes[1])
end

#=
# test if adding data works

# check dimensions
function _check_dimensions(vptdf::VirtualPTDF, v::Array{Float64})
    row_len = size(vptdf, 2)
    if length(array) != row_len
        error(
            "If using array as imput, it must be the same lenght as the number of column of Matrix BA",
        )
    end
end

# functions to set indeces
"""
Can allocate values (single number or entire row) in the PTDF matrix
(data[2]).
"""
function Base.setindex!(
    vptdf::VirtualPTDF,
    v::Union{Float64, Array{Float64}},
    row::Union{Integer, String},
    column::Union{Integer, Colon},
)

    # allocate line name
    if row isa String
        row = _get_line_index(vptdf, row)
    end
    vptdf.data[1][row] = vptdf.axes[1][row]

    # allocate value/row
    if length(v) > 1 && column isa Colon
        # check if dimensions match
        _check_dimensions(v)
        vptdf.axes[2][row] = v
    elseif length(v) == 1 && column isa Integer
        if isempty(vptdf.axes[2][row])
            error(
                "Single element of each row can be modified only if row is already present",
            )
        else
            vptdf.axes[2][row][column] = v
        end
    else
        error("value(s) to allocate must either be a non-empty Float64 or Array{FLoat64}")
    end
end

## tests #####################################################################

# set values

# first set values per each indexed line
ptdf_complete = PTDF(sys; linear_solver = "KLU")
ptdf_virtual = VirtualPTDF(sys)

# per each row
for i in axes(ptdf_complete.data, 1)
    for j in axes(ptdf_complete.data, 2)
        # check values
        setindex!(ptdf_virtual, ptdf_complete[i, j], i, j)
    end
end

# per each element
for i in axes(ptdf_complete.data, 1)
    for j in axes(ptdf_complete.data, 2)
        # check values
        setindex!(ptdf_virtual, ptdf_complete[i, j], i, j)
    end
end

# now check if they are the same
count_ = 0
for i in axes(ptdf_complete.data, 1)
    for j in axes(ptdf_complete.data, 2)
        # check values
        if !isapprox(ptdf_complete[i, j], ptdf_virtual[i, j]; atol = 1e-10)
            @show count_ = +1
        end
    end
end

# set values per each line name

# now check if they are the same

# get values
ptdf_complete = PTDF(sys; linear_solver = "KLU")
ptdf_virtual = VirtualPTDF(sys)

count_ = 0
for i in axes(ptdf_complete.data, 1)
    for j in axes(ptdf_complete.data, 2)
        # check values
        if !isapprox(ptdf_complete.data[i, j], ptdf_virtual[i, j]; atol = 1e-10)
            @show count_ = +1
            # @show ptdf_complete.data[i, j] - ptdf_virtual[i, j]
        end
    end
end

ptdf_virtual = VirtualPTDF(sys)
count1_ = 0
for i in axes(ptdf_complete.data, 1)
    # check values
    if !isapprox(ptdf_complete.data[i, :], ptdf_virtual[i, :]; atol = 1e-10)
        @show count1_ = +1
        # @show ptdf_complete.data[i, j] - ptdf_virtual[i, j]
    end
end

ptdf_virtual = VirtualPTDF(sys)
count2_ = 0
for (i, name) in enumerate(ptdf_complete.axes[1])
    if !isapprox(ptdf_complete.data[i, :], ptdf_virtual[name, :]; atol = 1e-10)
        @show count2_ = +1
        # @show ptdf_complete.data[i, j] - ptdf_virtual[i, j]
    end
end

"""
no differences
"""

# final check on entire matrix
new_ptdf = reduce(vcat, transpose.(ptdf_virtual.data[2]))
old_ptds = ptdf_complete.data
isapprox(new_ptdf, old_ptds; atol = 1e-10)
=#
