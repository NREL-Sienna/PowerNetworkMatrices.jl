"""
Power Transfer Distribution Factors (PTDF) indicate the incremental change in
real power that occurs on transmission lines due to real power injections
changes at the buses.

The PTDF struct is indexed using the Bus numbers and Branch names.

# Arguments
- `data<:AbstractArray{Float64, 2}`:
        the transposed PTDF matrix.
- `axes<:NTuple{2, Dict}`:
        Tuple containing two vectors: the first one showing the bus numbers,
        the second showing the branch names. The information contained in this
        field matches the axes of the fields `data`.
- `lookup<:NTuple{2, Dict}`:
        Tuple containing two dictionaries mapping the bus numbers and branch
        names with the indices of the matrix contained in `data`.
- `subnetworks::Dict{Int, Set{Int}}`:
        dictionary containing the set of bus indexes defining the subnetworks
        of the system.
- `tol::Base.RefValue{Float64}`:
        tolerance used for sparsifying the matrix (dropping element whose
        absolute value is below this threshold).
- `network_reduction::NetworkReduction`:
        Structure containing the details of the network reduction applied when computing the matrix
"""
struct PTDF{Ax, L <: NTuple{2, Dict}, M <: AbstractArray{Float64, 2}} <:
       PowerNetworkMatrix{Float64}
    data::M
    axes::Ax
    lookup::L
    subnetwork_axes::Dict{Int, Ax}
    tol::Base.RefValue{Float64}
    network_reduction_data::NetworkReductionData
end

get_axes(M::PTDF) = M.axes
get_lookup(M::PTDF) = M.lookup
get_ref_bus(M::PTDF) = sort!(collect(keys(M.subnetwork_axes)))
get_ref_bus_position(M::PTDF) = [get_bus_lookup(M)[x] for x in keys(M.subnetwork_axes)]
get_network_reduction_data(M::PTDF) = M.network_reduction_data
get_bus_axis(M::PTDF) = M.axes[1]
get_bus_lookup(M::PTDF) = M.lookup[1]
get_arc_axis(M::PTDF) = M.axes[2]
get_arc_lookup(M::PTDF) = M.lookup[2]

stores_transpose(::PTDF) = true

"""
Deserialize a PTDF from an HDF5 file.

# Arguments
- `filename::AbstractString`: File containing a serialized PTDF.
"""
PTDF(filename::AbstractString) = from_hdf5(PTDF, filename)

function _buildptdf_from_matrices(
    A::SparseArrays.SparseMatrixCSC{Int8, Int},
    BA::SparseArrays.SparseMatrixCSC{T, Int} where {T <: Union{Float32, Float64}},
    ref_bus_positions::Set{Int},
    dist_slack::Vector{Float64},
    linear_solver::String)
    if linear_solver == "KLU"
        PTDFm = _calculate_PTDF_matrix_KLU(A, BA, ref_bus_positions, dist_slack)
    elseif linear_solver == "Dense"
        # Convert SparseMatrices to Dense
        PTDFm = _calculate_PTDF_matrix_DENSE(
            A,
            BA,
            ref_bus_positions,
            dist_slack,
        )
    elseif linear_solver == "MKLPardiso"
        PTDFm =
            _calculate_PTDF_matrix_MKLPardiso(A, BA, ref_bus_positions, dist_slack)
    end

    return PTDFm
end

"""
Function for internal use only.

Computes the PTDF matrix by means of the KLU.LU factorization for sparse matrices.

# Arguments
- `A::SparseArrays.SparseMatrixCSC{Int8, Int}`:
        Incidence Matrix
- `BA::SparseArrays.SparseMatrixCSC{Float64, Int}`:
        BA matrix
- `ref_bus_positions::Set{Int}`:
        vector containing the indexes of the reference slack buses.
- `dist_slack::Vector{Float64}`:
        vector containing the weights for the distributed slacks.
"""
function _calculate_PTDF_matrix_KLU(
    A::SparseArrays.SparseMatrixCSC{Int8, Int},
    BA::SparseArrays.SparseMatrixCSC{Float64, Int},
    ref_bus_positions::Set{Int},
    dist_slack::Vector{Float64})
    linecount = size(BA, 2)
    buscount = size(BA, 1)

    ABA = calculate_ABA_matrix(A, BA, ref_bus_positions)
    K = klu(ABA)
    # initialize matrices for evaluation
    valid_ix = setdiff(1:buscount, ref_bus_positions)
    PTDFm_t = zeros(buscount, linecount)
    copyto!(PTDFm_t, BA)
    if !isempty(dist_slack) && length(ref_bus_positions) != 1
        error(
            "Distributed slack is not supported for systems with multiple reference buses.",
        )
    elseif isempty(dist_slack) && length(ref_bus_positions) < buscount
        PTDFm_t[valid_ix, :] = KLU.solve!(K, PTDFm_t[valid_ix, :])
        PTDFm_t[collect(ref_bus_positions), :] .= 0.0
        return PTDFm_t
    elseif length(dist_slack) == buscount
        @info "Distributed bus"
        PTDFm_t[valid_ix, :] = KLU.solve!(K, PTDFm_t[valid_ix, :])
        PTDFm_t[collect(ref_bus_positions), :] .= 0.0
        slack_array = dist_slack / sum(dist_slack)
        slack_array = reshape(slack_array, 1, buscount)
        return PTDFm_t .- (slack_array * PTDFm_t)
    else
        error("Distributed bus specification doesn't match the number of buses.")
    end

    return
end

function _binfo_check(binfo::Int)
    if binfo != 0
        if binfo < 0
            error("Illegal Argument in Inputs")
        elseif binfo > 0
            error("Singular value in factorization. Possibly there is an islanded bus")
        else
            @assert false
        end
    end
    return
end

"""
Function for internal use only.

Computes the PTDF matrix by means of the LAPACK and BLAS functions for dense matrices.

# Arguments
- `A::Matrix{Int8}`:
        Incidence Matrix
- `BA::Matrix{T} where {T <: Union{Float32, Float64}}`:
        BA matrix
- `ref_bus_positions::Set{Int}`:
        vector containing the indexes of the reference slack buses.
- `dist_slack::Vector{Float64})`:
        vector containing the weights for the distributed slacks.
"""
function _calculate_PTDF_matrix_DENSE(
    A::SparseArrays.SparseMatrixCSC{Int8, Int},
    BA::SparseArrays.SparseMatrixCSC{Float64, Int},
    ref_bus_positions::Set{Int},
    dist_slack::Vector{Float64})
    linecount = size(BA, 2)
    buscount = size(BA, 1)
    # Use dense calculation of ABA
    valid_ixs = setdiff(1:buscount, ref_bus_positions)
    ABA = Matrix(calculate_ABA_matrix(A, BA, ref_bus_positions))
    PTDFm_t = zeros(buscount, linecount)
    (ABA, bipiv, binfo) = getrf!(ABA)
    _binfo_check(binfo)
    BA = Matrix(BA[valid_ixs, :])
    if !isempty(dist_slack) && length(ref_bus_positions) != 1
        error(
            "Distributed slack is not supported for systems with multiple reference buses.",
        )
    elseif isempty(dist_slack) && length(ref_bus_positions) < buscount
        getrs!('N', ABA, bipiv, BA)
        PTDFm_t[valid_ixs, :] = BA
        return PTDFm_t
    elseif length(dist_slack) == buscount
        @info "Distributed bus"
        getrs!('N', ABA, bipiv, BA)
        PTDFm_t[valid_ixs, :] = BA
        slack_array = dist_slack / sum(dist_slack)
        slack_array = reshape(slack_array, 1, buscount)
        return PTDFm_t -
               gemm('N', 'N', ones(buscount, 1), gemm('N', 'N', slack_array, PTDFm_t))
    else
        error("Distributed bus specification doesn't match the number of buses.")
    end

    return
end

"""
Function for internal use only.

Computes the PTDF matrix by means of the MKL Pardiso for dense matrices.

# Arguments
- `A::SparseArrays.SparseMatrixCSC{Int8, Int}`:
        Incidence Matrix
- `BA::SparseArrays.SparseMatrixCSC{Float64, Int}`:
        BA matrix
- `ref_bus_positions::Set{Int}`:
        vector containing the indexes of the reference slack buses.
- `dist_slack::Vector{Float64}`:
        vector containing the weights for the distributed slacks.
"""
function _calculate_PTDF_matrix_MKLPardiso(
    A::SparseArrays.SparseMatrixCSC{Int8, Int},
    BA::SparseArrays.SparseMatrixCSC{Float64, Int},
    ref_bus_positions::Set{Int},
    dist_slack::Vector{Float64})
    linecount = size(BA, 2)
    buscount = size(BA, 1)

    ABA = calculate_ABA_matrix(A, BA, ref_bus_positions)
    @assert LinearAlgebra.issymmetric(ABA)
    ps = Pardiso.MKLPardisoSolver()
    Pardiso.set_matrixtype!(ps, Pardiso.REAL_SYM)
    Pardiso.pardisoinit(ps)
    # Pardiso.set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)
    defaults = Pardiso.get_iparms(ps)
    Pardiso.set_iparm!(ps, 1, 1)
    for (ix, v) in enumerate(defaults[2:end])
        Pardiso.set_iparm!(ps, ix + 1, v)
    end
    Pardiso.set_iparm!(ps, 2, 2)
    Pardiso.set_iparm!(ps, 59, 2)
    Pardiso.set_iparm!(ps, 6, 1)
    Pardiso.set_iparm!(ps, 12, 1)
    Pardiso.set_iparm!(ps, 11, 0)
    Pardiso.set_iparm!(ps, 13, 0)
    Pardiso.set_iparm!(ps, 32, 1)

    # initialize matrices for evaluation
    valid_ix = setdiff(1:buscount, ref_bus_positions)
    PTDFm_t = zeros(buscount, linecount)

    full_BA = Matrix(BA[valid_ix, :])
    if !isempty(dist_slack) && length(ref_bus_positions) != 1
        error(
            "Distributed slack is not supported for systems with multiple reference buses.",
        )
    elseif isempty(dist_slack) && length(ref_bus_positions) != buscount
        Pardiso.pardiso(ps, PTDFm_t[valid_ix, :], ABA, full_BA)
        PTDFm_t[valid_ix, :] = full_BA
        Pardiso.set_phase!(ps, Pardiso.RELEASE_ALL)
        Pardiso.pardiso(ps)
        return PTDFm_t
    elseif length(dist_slack) == buscount
        @info "Distributed bus"
        Pardiso.pardiso(ps, PTDFm_t[valid_ix, :], ABA, full_BA)
        PTDFm_t[valid_ix, :] = full_BA
        Pardiso.set_phase!(ps, Pardiso.RELEASE_ALL)
        Pardiso.pardiso(ps)
        slack_array = dist_slack / sum(dist_slack)
        slack_array = reshape(slack_array, 1, buscount)
        return PTDFm_t - ones(buscount, 1) * (slack_array * PTDFm_t)
    else
        error("Distributed bus specification doesn't match the number of buses.")
    end
    return
end

function PTDF(sys::PSY.System;
    dist_slack::Dict{Int, Float64} = Dict{Int, Float64}(),
    linear_solver = "KLU",
    tol::Float64 = eps(),
    network_reductions::Vector{NetworkReduction} = NetworkReduction[],
    kwargs...,
)
    Ymatrix = Ybus(
        sys;
        network_reductions = network_reductions,
        kwargs...,
    )
    A = IncidenceMatrix(Ymatrix)
    BA = BA_Matrix(Ymatrix)
    return PTDF(
        A,
        BA;
        dist_slack = dist_slack,
        linear_solver = linear_solver,
        tol = tol,
    )
end

"""
Builds the PTDF matrix from a system. The return is a PTDF array indexed with the bus numbers.

# Arguments
- `A::IncidenceMatrix`:
        Incidence Matrix (full structure)
- `BA::BA_Matrix`:
        BA matrix (full structure)

# Keyword Arguments
- `dist_slack::Vector{Float64}`:
        Vector of weights to be used as distributed slack bus.
        The distributed slack vector has to be the same length as the number of buses.
- `linear_solver::String`:
        Linear solver to be used. Options are "Dense", "KLU" and "MKLPardiso.
- `tol::Float64`:
        Tolerance to eliminate entries in the PTDF matrix (default eps()).
"""
function PTDF(
    A::IncidenceMatrix,
    BA::BA_Matrix;
    dist_slack::Dict{Int, Float64} = Dict{Int, Float64}(),
    linear_solver = "KLU",
    tol::Float64 = eps(),
)
    if !(isempty(dist_slack))
        dist_slack = redistribute_dist_slack(dist_slack, A, A.network_reduction_data)
    else
        dist_slack = Float64[]
    end
    validate_linear_solver(linear_solver)
    if !isequal(A.network_reduction_data, BA.network_reduction_data)
        error("A and BA matrices have non-equivalent network reductions.")
    end
    axes = BA.axes
    lookup = BA.lookup
    A_matrix = A.data
    subnetwork_axes = BA.subnetwork_axes
    ref_bus_positions = get_ref_bus_position(BA)
    S = _buildptdf_from_matrices(
        A_matrix,
        BA.data,
        Set(ref_bus_positions),
        dist_slack,
        linear_solver,
    )
    if tol > eps()
        return PTDF(
            sparsify(S, tol),
            axes,
            lookup,
            subnetwork_axes,
            Ref(tol),
            BA.network_reduction_data,
        )
    else
        return PTDF(
            S,
            axes,
            lookup,
            subnetwork_axes,
            Ref(tol),
            BA.network_reduction_data,
        )
    end
end

##############################################################################
########################### Auxiliary functions ##############################
##############################################################################

function Base.getindex(A::PTDF, branch_name::String, bus)
    multiplier, arc = get_branch_multiplier(A, branch_name)
    i, j = to_index(A, bus, arc)
    return A.data[i, j] * multiplier
end

# PTDF stores the transposed matrix. Overload indexing and how data is exported.
function Base.getindex(A::PTDF, arc, bus)
    i, j = to_index(A, bus, arc)
    return A.data[i, j]
end

function Base.getindex(
    A::PTDF,
    line_number::Union{Int, Colon},
    bus_number::Union{Int, Colon},
)
    return A.data[bus_number, line_number]
end

function get_ptdf_data(ptdf::PTDF)
    return transpose(ptdf.data)
end

function get_tol(ptdf::PTDF)
    return ptdf.tol
end

function redistribute_dist_slack(
    dist_slack::Dict{Int, Float64},
    A::IncidenceMatrix,
    nr::NetworkReductionData,
)
    dist_slack_vector = zeros(length(A.axes[2]))
    for (bus_no, dist_slack_factor) in dist_slack
        bus_no_ = get(nr.reverse_bus_search_map, bus_no, bus_no)
        if !haskey(A.lookup[2], bus_no_)
            throw(
                IS.InvalidValue(
                    "Bus number $bus_no_ not found in the incidence matrix. Correct your slack distribution specification.",
                ),
            )
        end
        ix = A.lookup[2][bus_no_]
        dist_slack_vector[ix] += dist_slack_factor
    end
    return dist_slack_vector
end
