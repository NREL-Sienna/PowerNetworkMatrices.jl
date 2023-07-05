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
struct VirtualLODF{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Float64}
    K::KLU.KLUFactorization{Float64, Int}
    BA::SparseArrays.SparseMatrixCSC{Float64, Int}
    A::SparseArrays.SparseMatrixCSC{Int8, Int}
    inv_PTDF_A_diag::Vector{Float64}
    ref_bus_positions::Set{Int}
    dist_slack::Vector{Float64}
    axes::Ax
    lookup::L
    valid_ix::Vector{Int}
    temp_data::Vector{Float64}
    cache::RowCache
    subnetworks::Dict{Int, Set{Int}}
    tol::Base.RefValue{Float64}
end

function _get_PTDF_A_diag(
    K::KLU.KLUFactorization{Float64, Int},
    BA::SparseArrays.SparseMatrixCSC{Float64, Int},
    A::SparseArrays.SparseMatrixCSC{Int8, Int64},
    ref_bus_positions::Set{Int},
    dist_slack::Vector{Float64},
)
    buscount, linecount = size(BA)
    # inizialize matrices for evaluation
    valid_ix = setdiff(1:buscount, ref_bus_positions)
    PTDFm_t = zeros(buscount, linecount)
    if !isempty(dist_slack) && length(ref_bus_positions) != 1
        error(
            "Distibuted slack is not supported for systems with multiple reference buses.",
        )
    elseif isempty(dist_slack) && length(ref_bus_positions) < buscount

        # get inverse of ABA
        Ix = SparseArrays.sparse(I, size(K, 1), size(K, 1))
        ABA_inv = zeros(Float64, size(Ix))
        ldiv!(ABA_inv, K, Ix)
        # multiply the matrix just for some elements
        diag_ = zeros(size(BA, 2))
        for i in 1:size(BA, 2) # per each column
            for j in BA.rowval[BA.colptr[i]:(BA.colptr[i + 1] - 1)]
                check_1 = sum(j .> ref_bus_positions)
                for k in BA.colptr[i]:(BA.colptr[i + 1] - 1)
                    if BA.rowval[k] ∉ ref_bus_positions && j ∉ ref_bus_positions
                        check_2 = sum(BA.rowval[k] .> ref_bus_positions)
                        diag_[i] +=
                            A[i, j] *
                            (BA.nzval[k] * ABA_inv[j - check_1, BA.rowval[k] - check_2])
                    end
                end
            end
        end

        return diag_

    elseif length(dist_slack) == buscount

        # TODO still to implement as above

        @info "Distributed bus"
        copyto!(PTDFm_t, BA)
        PTDFm_t[valid_ix, :] = KLU.solve!(K, PTDFm_t[valid_ix, :])
        PTDFm_t[ref_bus_positions, :] .= 0.0
        slack_array = dist_slack / sum(dist_slack)
        slack_array = reshape(slack_array, 1, buscount)
        return LinearAlgebra.diag(
            A * (PTDFm_t - (PTDFm_t * slack_array) * ones(1, buscount)),
        )
    else
        error("Distributed bus specification doesn't match the number of buses.")
    end
end

function VirtualLODF(
    sys::PSY.System; kwargs...,
)
    branches = get_ac_branches(sys)
    buses = get_buses(sys)

    return VirtualLODF(branches, buses; kwargs...)
end

function VirtualLODF(
    branches,
    buses::Vector{PSY.ACBus};
    dist_slack::Vector{Float64} = Float64[],
    tol::Float64 = eps(),
    max_cache_size::Int = MAX_CACHE_SIZE_MiB,
    persistent_lines::Vector{String} = String[],
)

    #Get axis names and lookups
    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in buses]
    axes = (line_ax, line_ax)
    # get matrices
    M, bus_ax_ref = calculate_adjacency(branches, buses)
    A, ref_bus_positions = calculate_A_matrix(branches, buses)
    BA = calculate_BA_matrix(branches, bus_ax_ref)
    K = klu(calculate_ABA_matrix(A, BA, ref_bus_positions))
    # get lookups, reference bus positions and subnetworks
    line_ax_ref = make_ax_ref(line_ax)
    look_up = (line_ax_ref, line_ax_ref)
    ref_bus_positions = find_slack_positions(buses)
    subnetworks = find_subnetworks(M, bus_ax)
    # check subnetworks
    if length(subnetworks) > 1
        @info "Network is not connected, using subnetworks"
        subnetworks = assing_reference_buses(subnetworks, ref_bus_positions)
    end
    # get diagonal of PTDF
    PTDF_diag = _get_PTDF_A_diag(
        K,
        BA,
        A,
        ref_bus_positions,
        dist_slack,
    )
    PTDF_diag[PTDF_diag .> 1 - 1e-6] .= 0.0
    # initialize structure
    temp_data = zeros(length(bus_ax))
    valid_ix = setdiff(1:length(bus_ax), ref_bus_positions)
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

    return VirtualLODF(
        K,
        BA,
        A,
        1.0 ./ (1.0 .- PTDF_diag),
        ref_bus_positions,
        dist_slack,
        axes,
        look_up,
        valid_ix,
        temp_data,
        empty_cache,
        subnetworks,
        Ref(tol),
    )
end

# Overload Base functions

"""
Checks if the any of the fields of VirtualPTDF is empty.
"""
# ! does not check as expected: if first line is true it does not check the others
function Base.isempty(vlodf::VirtualLODF)
    !isempty(vlodf.K.L) && return false
    !isempty(vlodf.K.U) && return false
    !isempty(vlodf.BA) && return false
    !isempty(vlodf.A) && return false
    !isempty(vlodf.PTDF_A_diag) && return false
    !isempty(vlodf.ref_bus_positions) && return false
    !isempty(vlodf.dist_slack) && return false
    !isempty(vlodf.axes) && return false
    !isempty(vlodf.lookup) && return false
    !isempty(vlodf.valid_ix) && return false
    !isempty(vlodf.cache) && return false
    !isempty(vlodf.temp_data) && return false
    !isempty(subnetworks) && return false
    !isempty(tol) && return false
    return true
end

"""
Gives the size of the whole PTDF matrix, not the number of rows stored.
"""
# ! test
Base.size(vlodf::VirtualLODF) = (size(vlodf.BA, 2), size(vlodf.BA, 2))
"""
Gives the cartesian indexes of the PTDF matrix (same as the BA one).
"""
Base.eachindex(vlodf::VirtualLODF) = CartesianIndices(size(vlodf))

if isdefined(Base, :print_array) # 0.7 and later
    Base.print_array(io::IO, X::VirtualLODF) = "VirtualLODF"
end

function _getindex(
    vlodf::VirtualLODF,
    row::Int,
    column::Union{Int, Colon},
)
    # check if value is in the cache
    if haskey(vlodf.cache, row)
        return vlodf.cache[row][column]
    else
        # evaluate the value for the PTDF column
        # Needs improvement (not much found...)

        lin_solve = KLU.solve!(vlodf.K, Vector(vlodf.BA[vlodf.valid_ix, row]))
        # get full ptdf row
        for i in eachindex(vlodf.valid_ix)
            vlodf.temp_data[vlodf.valid_ix[i]] = lin_solve[i]
        end
        lodf_row = (vlodf.A * vlodf.temp_data) .* vlodf.inv_PTDF_A_diag
        lodf_row[row] = -1.0
        # add slack bus value (zero) and make copy of temp into the cache
        if get_tol(vlodf) > eps()
            vlodf.cache[row] = make_entries_zero!(deepcopy(lodf_row), get_tol(vlodf))
        else
            vlodf.cache[row] = deepcopy(lodf_row)
        end
        return vlodf.cache[row][column]
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
function Base.getindex(vlodf::VirtualLODF, row, column)
    row_, column_ = to_index(vlodf, row, column)
    return _getindex(vlodf, row_, column_)
end

# Define for ambiguity resolution
function Base.getindex(vlodf::VirtualLODF, row::Integer, column::Integer)
    return _getindex(vlodf, row, column)
end

"""
!!! STILL TO IMPLEMENT !!!
"""
Base.setindex!(::VirtualLODF, _, idx...) = error("Operation not supported by VirtualPTDF")

"""
!!! STILL TO IMPLEMENT !!!
"""
Base.setindex!(::VirtualLODF, _, ::CartesianIndex) =
    error("Operation not supported by VirtualPTDF")

"""
PTDF data is stored in the the cache
it is a nested vector containing an array for the names of each row,
the PTDF's matrices rows and how many times they were evaluated
"""

# ! change it so to get only the non-empty values

get_data(mat::VirtualLODF) = mat.cache

function get_branch_ax(ptdf::VirtualLODF)
    return ptdf.axes[1]
end

""" Gets the tolerance used for sparsifying the rows of the PTDF matrix"""
function get_tol(vptdf::Union{VirtualPTDF, VirtualLODF})
    return vptdf.tol[]
end
