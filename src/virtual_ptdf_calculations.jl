"""
Power Transfer Distribution Factors (PTDF) indicate the incremental change in
real power that occurs on transmission lines due to real power injections
changes at the buses.

The PTDF struct is indexed using the Bus numbers and branch names.

# Arguments
- `K::KLU.KLUFactorization{Float64, Int}`:
        LU factorization matrices of the ABA matrix, evaluated by means of KLU
- `BA::SparseArrays.SparseMatrixCSC{Float64, Int}`:
        BA matric
- `ref_bus_positions::Set{Int}`:
        Vector containing the indexes of the columns of the BA matrix corresponding
        to the refence buses
- `dist_slack::Vector{Float64}`:
        Vector of weights to be used as distributed slack bus.
        The distributed slack vector has to be the same length as the number of buses.
- `axes<:NTuple{2, Dict}`:
        Tuple containing two vectors (the first one showing the branches names,
        the second showing the buses numbers).
- `lookup<:NTuple{2, Dict}`:
        Tuple containing two dictionaries, the first mapping the branches
        and buses with their enumerated indexes.
- `temp_data::Vector{Float64}`:
        temporary vector for internal use.
- `cache::RowCache`:
        cache were PTDF rows are stored.
- `subnetworks::Dict{Int, Set{Int}}`:
        dictionary containing the subsets of buses defining the different subnetwork of the system.
- `tol::Base.RefValue{Float64}`:
        tolerance related to scarification and values to drop.
"""
struct VirtualPTDF{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Float64}
    K::KLU.KLUFactorization{Float64, Int}
    BA::SparseArrays.SparseMatrixCSC{Float64, Int}
    ref_bus_positions::Set{Int}
    dist_slack::Vector{Float64}
    axes::Ax
    lookup::L
    temp_data::Vector{Float64}
    valid_ix::Vector{Int}
    cache::RowCache
    subnetworks::Dict{Int, Set{Int}}
    tol::Base.RefValue{Float64}
end

"""
Builds the PTDF matrix from a group of branches and buses. The return is a
PTDF array indexed with the branch numbers.

# Arguments
- `branches`:
        Vector of the system's AC branches.
- `buses::Vector{PSY.Bus}`:
        Vector of the system's buses.
- `dist_slack::Vector{Float64} = Float64[]`:
        Vector of weights to be used as distributed slack bus.
        The distributed slack vector has to be the same length as the number of buses.
- `tol::Float64 = eps()`:
        Tolerance related to sparsification and values to drop.
- `max_cache_size::Int`:
        max cache size in MiB (inizialized as MAX_CACHE_SIZE_MiB).
- `persistent_lines::Vector{String}`:
        line to be evaluated as soon as the VirtualPTDF is created (initialized as empty vector of strings).
"""
function VirtualPTDF(
    branches,
    buses::Vector{PSY.Bus};
    dist_slack::Vector{Float64} = Float64[],
    tol::Float64 = eps(),
    max_cache_size::Int = MAX_CACHE_SIZE_MiB,
    persistent_lines::Vector{String} = String[])
    if length(dist_slack) != 0
        @info "Distributed bus"
    end

    #Get axis names
    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in buses]
    axes = (line_ax, bus_ax)
    M, bus_ax_ref = calculate_adjacency(branches, buses)
    line_ax_ref = make_ax_ref(line_ax)
    look_up = (line_ax_ref, bus_ax_ref)
    A, ref_bus_positions = calculate_A_matrix(branches, buses)
    BA = calculate_BA_matrix(branches, bus_ax_ref)
    ABA = calculate_ABA_matrix(A, BA, ref_bus_positions)
    ref_bus_positions = find_slack_positions(buses)
    subnetworks = find_subnetworks(M, bus_ax)
    if length(subnetworks) > 1
        @info "Network is not connected, using subnetworks"
        subnetworks = assing_reference_buses(subnetworks, ref_bus_positions)
    end
    temp_data = zeros(length(bus_ax))
    if isempty(persistent_lines)
        empty_cache =
            RowCache(max_cache_size * MiB, Set{Int}(), length(bus_ax) * sizeof(Float64))
    else
        init_persistent_dict = Set{Int}(line_ax_ref[k] for k in persistent_lines)
        empty_cache =
            RowCache(
                max_cache_size * MiB,
                init_persistent_dict,
                length(bus_ax) * sizeof(Float64),
            )
    end
    return VirtualPTDF(
        klu(ABA),
        BA,
        ref_bus_positions,
        dist_slack,
        axes,
        look_up,
        temp_data,
        setdiff(1:length(temp_data), ref_bus_positions),
        empty_cache,
        subnetworks,
        Ref(tol),
    )
end

"""
Builds the Virtual PTDF matrix from a system. The return is a VirtualPTDF
struct with an empty cache.
"""
function VirtualPTDF(
    sys::PSY.System; kwargs...)
    branches = get_ac_branches(sys)
    buses = get_buses(sys)

    return VirtualPTDF(branches, buses; kwargs...)
end

# Overload Base functions

"""
Checks if the any of the fields of VirtualPTDF is empty.
"""
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

"""
Gives the size of the whole PTDF matrix, not the number of rows stored.
"""
Base.size(vptdf::VirtualPTDF) = size(vptdf.BA)
"""
Gives the cartesian indexes of the PTDF matrix (same as the BA one).
"""
Base.eachindex(vptdf::VirtualPTDF) = CartesianIndices(size(vptdf.BA))

if isdefined(Base, :print_array) # 0.7 and later
    Base.print_array(io::IO, X::VirtualPTDF) = "VirtualPTDF"
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
        # Needs improvement
        valid_ix = vptdf.valid_ix
        lin_solve = KLU.solve!(vptdf.K, Vector(vptdf.BA[valid_ix, row]))
        buscount = size(vptdf, 1)

        if !isempty(vptdf.dist_slack) && length(vptdf.ref_bus_positions) != 1
            error(
                "Distibuted slack is not supported for systems with multiple reference buses.",
            )
        elseif isempty(vptdf.dist_slack) && length(vptdf.ref_bus_positions) < buscount
            for i in eachindex(valid_ix)
                vptdf.temp_data[valid_ix[i]] = lin_solve[i]
            end
            vptdf.cache[row] = deepcopy(vptdf.temp_data)
        elseif length(vptdf.dist_slack) == buscount
            for i in eachindex(valid_ix)
                vptdf.temp_data[valid_ix[i]] = lin_solve[i]
            end
            slack_array = vptdf.dist_slack / sum(vptdf.dist_slack)
            slack_array = reshape(slack_array, buscount)
            vptdf.cache[row] =
                deepcopy(vptdf.temp_data .- dot(vptdf.temp_data, slack_array))
        else
            error("Distributed bus specification doesn't match the number of buses.")
        end

        # add slack bus value (zero) and make copy of temp into the cache
        if get_tol(vptdf) > eps()
            # vptdf.cache[row] = make_entries_zero!(deepcopy(vptdf.temp_data), get_tol(vptdf))
            make_entries_zero!(vptdf.cache[row], get_tol(vptdf))
        end

        return vptdf.cache[row][column]
    end
end

"""
Gets the value of the element of the PTDF matrix given the row and column indices
corresponding to the branch and buses one respectively. If `column` is a Colon then
the entire row is returned.

# Arguments
- `vptdf::VirtualPTDF`:
        VirtualPTDF struct where to evaluate and store the values.
- `row`:
        Branch index.
- `column`:
        Bus index. If Colon then get the values of the whole row.
"""
function Base.getindex(vptdf::VirtualPTDF, row, column)
    row_, column_ = to_index(vptdf, row, column)
    return _getindex(vptdf, row_, column_)
end

# Define for ambiguity resolution
function Base.getindex(vptdf::VirtualPTDF, row::Integer, column::Integer)
    return _getindex(vptdf, row, column)
end

"""
!!! STILL TO IMPLEMENT !!!
"""
Base.setindex!(::VirtualPTDF, _, idx...) = error("Operation not supported by VirtualPTDF")

"""
!!! STILL TO IMPLEMENT !!!
"""
Base.setindex!(::VirtualPTDF, _, ::CartesianIndex) =
    error("Operation not supported by VirtualPTDF")

"""
PTDF data is stored in the the cache
it is a nested vector containing an array for the names of each row,
the PTDF's matrices rows and how many times they were evaluated
"""

# ! change it so to get only the non-empty values

get_data(mat::VirtualPTDF) = mat.cache

function get_branch_ax(ptdf::VirtualPTDF)
    return ptdf.axes[1]
end

function get_bus_ax(ptdf::VirtualPTDF)
    return ptdf.axes[2]
end
