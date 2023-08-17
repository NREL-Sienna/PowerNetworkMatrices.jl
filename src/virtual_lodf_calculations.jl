"""
The Virtual Line Outage Distribution Factor (VirtualLODF) structure gathers
the rows of the LODF matrix as they are evaluated on-the-go. These rows are 
evalauted independently, cached in the structure and do not require the 
computation of the whole matrix (therefore significantly reducing the 
computational requirements).

The VirtualLODF is initialized with no row stored.

The VirtualLODF struct is indexed using branch names.

# Arguments
- `K::KLU.KLUFactorization{Float64, Int}`:
        LU factorization matrices of the ABA matrix, evaluated by means of KLU.
- `BA::SparseArrays.SparseMatrixCSC{Float64, Int}`:
        BA matrix.
- `A::SparseArrays.SparseMatrixCSC{Int8, Int}`:
        Incidence matrix.
- `inv_PTDF_A_diag::Vector{Float64}`:
        Vector contiaining the element-wise reciprocal of the diagonal elements
        coming from multuiplying the PTDF matrix with th Incidence matrix
- `ref_bus_positions::Set{Int}`:
        Vector containing the indexes of the rows of the transposed BA matrix 
        corresponding to the refence buses.
- `axes<:NTuple{2, Dict}`:
        Tuple containing two vectors showing the branch names.
- `lookup<:NTuple{2, Dict}`:
        Tuple containing two dictionaries, mapping the branches names 
        the enumerated row indexes indexes.
- `valid_ix::Vector{Int}`:
        Vector containing the row/columns indices of matrices related the buses
        which are not slack ones.
- `temp_data::Vector{Float64}`:
        Temporary vector for internal use.
- `cache::RowCache`:
        Cache were LODF rows are stored.
- `subnetworks::Dict{Int, Set{Int}}`:
        Dictionary containing the subsets of buses defining the different subnetwork of the system.
- `tol::Base.RefValue{Float64}`:
        Tolerance related to scarification and values to drop.
"""
struct VirtualLODF{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Float64}
    K::KLU.KLUFactorization{Float64, Int}
    BA::SparseArrays.SparseMatrixCSC{Float64, Int}
    A::SparseArrays.SparseMatrixCSC{Int8, Int}
    inv_PTDF_A_diag::Vector{Float64}
    ref_bus_positions::Set{Int}
    axes::Ax
    lookup::L
    valid_ix::Vector{Int}
    temp_data::Vector{Float64}
    cache::RowCache
    subnetworks::Dict{Int, Set{Int}}
    tol::Base.RefValue{Float64}
end

"""
Sets to zero those elements of each LODF matrix row whose absolute values are 
below the threshold specified by the field "tol".

# Arguments
- `mat::VirtualLODF`:
        VirtualLODF structure
- `tol::Float64`:
        tolerance
"""
function drop_small_entries!(mat::VirtualLODF, tol::Float64)
    if tol < mat.tol[]
        @info "Specified tolerance is smaller than the current tolerance."
    end
    for i in keys(mat.cache.temp_cache)
        make_entries_zero!(mat.cache[i], tol)
    end
    mat.tol[] = tol
    return
end

function _get_PTDF_A_diag(
    K::KLU.KLUFactorization{Float64, Int},
    BA::SparseArrays.SparseMatrixCSC{Float64, Int},
    A::SparseArrays.SparseMatrixCSC{Int8, Int64},
    ref_bus_positions::Set{Int},
)
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
    buses::Vector{PSY.Bus};
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
    temp_data = zeros(length(bus_ax))
    valid_ix = setdiff(1:length(bus_ax), ref_bus_positions)
    PTDF_diag = _get_PTDF_A_diag(
        K,
        BA,
        A,
        ref_bus_positions,
    )
    PTDF_diag[PTDF_diag .> 1 - 1e-6] .= 0.0
    # initialize structure
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
Checks if the any of the fields of VirtualLODF is empty.
"""
function Base.isempty(vlodf::VirtualLODF)
    !isempty(vlodf.K.L) && return false
    !isempty(vlodf.K.U) && return false
    !isempty(vlodf.BA) && return false
    !isempty(vlodf.A) && return false
    !isempty(vlodf.inv_PTDF_A_diag) && return false
    !isempty(vlodf.ref_bus_positions) && return false
    !isempty(vlodf.axes) && return false
    !isempty(vlodf.lookup) && return false
    !isempty(vlodf.cache) && return false
    !isempty(subnetworks) && return false
    !isempty(tol) && return false
    return true
end

"""
Shows the size of the whole LODF matrix, not the number of rows stored.
"""
# ! test
Base.size(vlodf::VirtualLODF) = (size(vlodf.BA, 2), size(vlodf.BA, 2))
"""
Gives the cartesian indexes of the LODF matrix.
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

        # evaluate the value for the LODF column
        
        # TODO: needs improvement to speed up computation (not much found...)

        lin_solve = KLU.solve!(vlodf.K, Vector(vlodf.BA[vlodf.valid_ix, row]))
        # get full lodf row
        for i in eachindex(vlodf.valid_ix)
            vlodf.temp_data[vlodf.valid_ix[i]] = lin_solve[i]
        end

        # now get the LODF row
        lodf_row = (vlodf.A * vlodf.temp_data) .* vlodf.inv_PTDF_A_diag
        lodf_row[row] = -1.0

        if get_tol(vlodf) > eps()
            vlodf.cache[row] = deepcopy(sparsify(lodf_row, get_tol(vlodf)))
        else
            vlodf.cache[row] = deepcopy(lodf_row)
        end
        return vlodf.cache[row][column]
    end
end

"""
Gets the value of the element of the LODF matrix given the row and column indices
corresponding to the selected and outage branch respectively. If `column` is a Colon then
the entire row is returned.

# Arguments
- `vlodf::VirtualLODF`:
        VirtualLODF struct where to evaluate and store the row values.
- `row`:
        selected line name
- `column`:
        outage line name. If `Colon` then get the values of the whole row.
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
Base.setindex!(::VirtualLODF, _, idx...) = error("Operation not supported by VirtualLODF")

"""
!!! STILL TO IMPLEMENT !!!
"""
Base.setindex!(::VirtualLODF, _, ::CartesianIndex) =
    error("Operation not supported by VirtualLODF")

get_lodf_data(mat::VirtualLODF) = mat.cache.temp_cache

function get_branch_ax(mat::VirtualLODF)
    return mat.axes[1]
end

""" Gets the tolerance used for sparsifying the rows of the VirtualLODF matrix"""
function get_tol(mat::VirtualLODF)
    return mat.tol[]
end
