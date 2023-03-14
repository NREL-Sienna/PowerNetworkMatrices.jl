"""
Power Transfer Distribution Factors (PTDF) indicate the incremental change in real power that occurs on transmission lines due to real power injections changes at the buses.

The PTDF struct is indexed using the Bus numbers and branch names
"""
struct VirtualPTDF{Ax, L <: NTuple{2, Dict}, T <: Real} <: PowerNetworkMatrix{Real}
    K::KLU.KLUFactorization{Float64, Int32}
    BA::SparseArrays.SparseMatrixCSC{T, Int32}
    slack_positions::Vector{Int64}
    dist_slack::Vector{Float64}
    axes::Ax
    lookup::L
    data::Dict{String, Array{Float64}}
    tol::Float64
end

function _build_virtualptdf(
    branches,
    nodes::Vector{PSY.Bus})
    slack_positions = find_slack_positions(nodes)
    A, _ = calculate_A_matrix(branches, nodes)
    BA = calculate_BA_matrix(branches, slack_positions, _make_ax_ref(nodes))
    ABA = calculate_ABA_matrix(A, BA, slack_positions)
    K = klu(ABA)

    return BA, K , A, slack_positions
end

function _buildptdf_row(
    K::KLU.KLUFactorization{Float64, Int32},
    BA_row::SparseArrays.SparseVector{Float64, Int32})
    # inizialize matrices for evaluation
    col = zeros(length(BA_row))
    ldiv!(col, K, BA_row)

    return col
end

# !
"""
Builds the PTDF matrix from a group of branches and nodes. The return is a PTDF array indexed with the bus numbers.

# Keyword arguments
- `dist_slack::Vector{Float64}`: Vector of weights to be used as distributed slack bus.
    The distributed slack vector has to be the same length as the number of buses
"""
# !
function VirtualPTDF(
    branches,
    nodes::Vector{PSY.Bus},
    dist_slack::Vector{Float64} = [0.1],
    tol::Float64 = 1e-7)

    #Get axis names
    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in nodes]
    axes = (line_ax, bus_ax)
    look_up = (_make_ax_ref(line_ax), _make_ax_ref(bus_ax))
    BA, K, _, slack_positions = _build_virtualptdf(branches, nodes)
    empty_cache = Dict{String, Array{Flaot64}}()
    return VirtualPTDF(K, BA, slack_positions, dist_slack, axes, look_up, empty_cache, tol)

end

# !
"""
Builds the PTDF matrix from a system. The return is a PTDF array indexed with the bus numbers.

# Keyword arguments
- `dist_slack::Vector{Float64}`: Vector of weights to be used as distributed slack bus.
    The distributed slack vector has to be the same length as the number of buses
"""
# !
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
    !isempty(vptdf.data) && return false
    return true
end
# Size related overload will work on the BA matrix
Base.size(vptdf::VirtualPTDF) = size(vptdf.BA)
Base.eachindex(vptdf::VirtualPTDF) = CartesianIndices(size(vptdf.BA))

# Get indices
function _get_line_index(vptdf::VirtualPTDF, row::String)
    row_ = findall(x->x==row, vptdf.axes[1])
    if length(row_) > 1
        error("multiple lines with the same name $row in vptdf.axes[1]")
    elseif length(row_) > 1
        error("no line with name $row in vptdf.axes[1]")
    else
        row = row_[1]
    end
    return row
end

function _get_value(vptdf::VirtualPTDF, branch::String, column::Union{Int, Colon})
    # get stored value if present in the dictionarr
    if branch in keys(vptdf.data)
        return vptdf.data[2][row]
    end
end

function Base.getindex(vptdf::VirtualPTDF, row::Union{Int, String}, column::Union{Int, Colon})

    # ! add tolerance field --> initialize as Float64 (round number below this to zero)
    # ! it is imported by vptdf

    # at first get the index value if needed
    if row isa String
        row = _get_line_index(vptdf, row)
    end
    # check if value is in the cache
    value_ = _get_value(vptdf, row, column)
    if !isempty(value_)
        return value_
    else
        # evaluate the value for the PTDF column
        PTDF_row = _buildptdf_row(vptdf.K, vptdf.BA[row, :])
        # add slack bus value (zero)
        PTDF_row = vcat(PTDF_row[1:vptdf.slack_positions[1]-1], [0], PTDF_row[vptdf.slack_positions[1]:end])
        # add the valute in cache
        vptdf.data[1][row] = vptdf.axes[1][row]
        vptdf.data[2][row] = PTDF_row
        # return value
        if column isa Colon
            return PTDF_row
        else
            return PTDF_row[column]
        end
    end
end

# Get indices
# Base.setindex!(A::VirtualPTDF, v, idx...) = A.data[to_index(A, idx...)...] = v
# Base.setindex!(A::VirtualPTDF, v, idx::CartesianIndex) = A.data[idx] = v

""" PTDF data is stored in the the cache
    it is a nested vecotr containing an array for the names of each row,
    the PTDF's matrices rows and how many times they were evaluated """

# ! change it so to get only the non-empty values

get_data(mat::PowerNetworkMatrix) = mat.data