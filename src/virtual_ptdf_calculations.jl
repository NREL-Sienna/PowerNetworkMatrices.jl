# TODO: add descriptions to functions

"""
Power Transfer Distribution Factors (PTDF) indicate the incremental change in real power that occurs on transmission lines due to real power injections changes at the buses.

The PTDF struct is indexed using the Bus numbers and branch names
"""
struct VirtualPTDF{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Float64}
    K::KLU.KLUFactorization{Float64, Int}
    BA::SparseArrays.SparseMatrixCSC{Float64, Int}
    ref_bus_positions::Vector{Int}
    dist_slack::Vector{Float64}
    axes::Ax
    lookup::L
    temp_data::Vector{Float64}
    cache::Dict{Int, Array{Float64}}
    tol::Base.RefValue{Float64}
end

"""
Builds the PTDF matrix from a group of branches and nodes. The return is a PTDF array indexed with the bus numbers.

# Keyword arguments
- `dist_slack::Vector{Float64}`: Vector of weights to be used as distributed slack bus.
    The distributed slack vector has to be the same length as the number of buses
"""
function VirtualPTDF(
    branches,
    nodes::Vector{PSY.Bus},
    dist_slack::Vector{Float64} = [0.1],
    tol::Float64 = eps())

    #Get axis names
    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in nodes]
    axes = (line_ax, bus_ax)
    look_up = (make_ax_ref(line_ax), make_ax_ref(bus_ax))
    A, ref_bus_positions = calculate_A_matrix(branches, nodes)
    BA = calculate_BA_matrix(branches, ref_bus_positions, make_ax_ref(nodes))
    ABA = calculate_ABA_matrix(A, BA, ref_bus_positions)
    empty_cache = Dict{Int, Array{Float64}}()
    temp_data = zeros(length(bus_ax))
    return VirtualPTDF(
        klu(ABA),
        BA,
        ref_bus_positions,
        dist_slack,
        axes,
        look_up,
        temp_data,
        empty_cache,
        Ref(tol),
    )
end

"""
Builds the PTDF matrix from a system. The return is a PTDF array indexed with the bus numbers.

# Keyword arguments
- `dist_slack::Vector{Float64}`: Vector of weights to be used as distributed slack bus.
    The distributed slack vector has to be the same length as the number of buses
"""
function VirtualPTDF(
    sys::PSY.System,
    dist_slack::Vector{Float64} = [0.1])
    branches = get_ac_branches(sys)
    nodes = get_buses(sys)

    return VirtualPTDF(branches, nodes, dist_slack)
end

# Overload Base functions
function Base.isempty(vptdf::VirtualPTDF)
    !isempty(vptdf.K.L) && return false
    !isempty(vptdf.K.U) && return false
    !isempty(vptdf.BA) && return false
    !isempty(vptdf.dist_slack) && return false
    !isempty(vptdf.axes) && return false
    !isempty(vptdf.lookup) && return false
    !isempty(vptdf.cache) && return false
    !isempty(vptdf.temp_data) && return false
    return true
end
# Size related overload will work on the BA matrix
Base.size(vptdf::VirtualPTDF) = size(vptdf.BA)
Base.eachindex(vptdf::VirtualPTDF) = CartesianIndices(size(vptdf.BA))

if isdefined(Base, :print_array) # 0.7 and later
    Base.print_array(io::IO, X::VirtualPTDF) = "VirtualPTDF"
end

# Get indices
function _get_line_index(vptdf::VirtualPTDF, row::String)
    row_ = findall(x -> x == row, vptdf.axes[1])
    if length(row_) > 1
        error("multiple lines with the same name $row in vptdf.axes[1]")
    elseif length(row_) > 1
        error("no line with name $row in vptdf.axes[1]")
    else
        row = row_[1]
    end
    return row
end

function _get_value(vptdf::VirtualPTDF, row::Int)
    # get stored value if present in the dictionarr
    if branch in keys(vptdf.data)
        return vptdf.data[2][row]
    end
end

function Base.getindex(vptdf::VirtualPTDF, row, column)
    row_, column_ = to_index(vptdf, row, column)
    return _getindex(vptdf, row_, column_)
end

# Define for ambiguity resolution
function Base.getindex(vptdf::VirtualPTDF, row::Integer, column::Integer)
    return _getindex(vptdf, row, column)
end

function _getindex(
    vptdf::VirtualPTDF,
    row::Int,
    column::Union{Int, Colon},
)

    # check if value is in the cache
    if haskey(vptdf.cache, row)
        return vptdf.cache[row][column]
    else
        # evaluate the value for the PTDF column
        vptdf.temp_data[setdiff(1:end, vptdf.ref_bus_positions)] .= KLU.solve!(vptdf.K, Vector(vptdf.BA[row, :]))
        # add slack bus value (zero) and make copy of temp into the cache
        vptdf.cache[row] = deepcopy(vptdf.temp_data)
        return vptdf.cache[row][column]
    end
end

Base.setindex!(::VirtualPTDF, _, idx...) = error("Operation not supported by VirtualPTDF")
Base.setindex!(::VirtualPTDF, _, ::CartesianIndex) =
    error("Operation not supported by VirtualPTDF")

""" PTDF data is stored in the the cache
    it is a nested vecotr containing an array for the names of each row,
    the PTDF's matrices rows and how many times they were evaluated """

# ! change it so to get only the non-empty values

get_data(mat::VirtualPTDF) = mat.cache
