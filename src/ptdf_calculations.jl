"""
Power Transfer Distribution Factors (PTDF) indicate the incremental change in real power that occurs on transmission lines due to real power injections changes at the buses.

The PTDF struct is indexed using the Bus numbers and Branch names.

# Arguments
- `data<:AbstractArray{Float64, 2}`:
        the actual Incidence matrix.
- `axes<:NTuple{2, Dict}`:
        Tuple containing two vectors (the first one showing the branches names,
        the second showing the buses numbers).
- `lookup<:NTuple{2, Dict}`:
        Tuple containing two dictionaries, the first mapping the branches
        and buses with their enumerated indexes.
- `subnetworks::Dict{Int, Set{Int}}`:
        dictionary containing the set of bus indexes defining the subnetworks
        of the system.
- `tol::Base.RefValue{Float64}`:
        tolerance used for sportifying the matrix (dropping element whose
        absolute value is below this threshold).
"""
struct PTDF{Ax, L <: NTuple{2, Dict}, M <: AbstractArray{Float64, 2}} <:
       PowerNetworkMatrix{Float64}
    data::M
    axes::Ax
    lookup::L
    subnetworks::Dict{Int, Set{Int}}
    ref_bus_positions::Set{Int}
    tol::Base.RefValue{Float64}
end

"""
Deserialize a PTDF from an HDF5 file.

# Arguments
- `filename::AbstractString`: File containing a serialized PTDF.
"""
PTDF(filename::AbstractString) = from_hdf5(PTDF, filename)

"""
Elements whose values are below "tol" are set to zero (dropped in case of Sparse matrix).
"""
function drop_small_entries!(mat::PTDF, tol::Float64)
    if tol > mat.tol[]
        @info "Specified tolerance is smaller than the current tolerance."
    end
    make_entries_zero!(mat.data, tol)
    mat.tol[] = tol
    return
end

"""
Takes and existing PTDF and makes it sparse given a certain tolerance "tol".
"""
function make_sparse_PTDF(mat::PTDF{Ax, L, Matrix{Float64}}, tol::Float64) where {Ax, L}
    new_mat = sparsify(mat.data, tol)
    return PTDF(new_mat, mat.axes, mat.lookup, ref_bus_positions, mat.subnetworks, Ref(tol))
end

function _buildptdf(
    branches,
    nodes::Vector{PSY.Bus},
    bus_lookup::Dict{Int, Int},
    dist_slack::Vector{Float64},
    linear_solver::String)
    if linear_solver == "KLU"
        PTDFm, A = calculate_PTDF_matrix_KLU(
            branches,
            nodes,
            bus_lookup,
            dist_slack,
        )
    elseif linear_solver == "Dense"
        PTDFm, A = calculate_PTDF_matrix_DENSE(
            branches,
            nodes,
            bus_lookup,
            dist_slack,
        )
    elseif linear_solver == "MKLPardiso"
        PTDFm, A = calculate_PTDF_matrix_MKLPardiso(
            branches,
            nodes,
            bus_lookup,
            dist_slack,
        )
    end

    return PTDFm, A
end

function _buildptdf_from_matrices(
    A::IncidenceMatrix,
    BA::SparseArrays.SparseMatrixCSC{T, Int} where {T <: Union{Float32, Float64}},
    dist_slack::Vector{Float64},
    linear_solver::String)
    if linear_solver == "KLU"
        PTDFm = _calculate_PTDF_matrix_KLU(A.data, BA, A.ref_bus_positions, dist_slack)
    elseif linear_solver == "Dense"
        # Convert SparseMatrices to Dense
        PTDFm = _calculate_PTDF_matrix_DENSE(
            Matrix(A.data),
            Matrix(BA),
            A.ref_bus_positions,
            dist_slack,
        )
    elseif linear_solver == "MKLPardiso"
        PTDFm =
            _calculate_PTDF_matrix_MKLPardiso(A.data, BA, A.ref_bus_positions, dist_slack)
    end

    return PTDFm
end

"""
Funciton for internal use only.

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
    # inizialize matrices for evaluation
    valid_ix = setdiff(1:buscount, ref_bus_positions)
    PTDFm_t = zeros(buscount, linecount)
    if !isempty(dist_slack) && length(ref_bus_positions) != 1
        error(
            "Distibuted slack is not supported for systems with multiple reference buses.",
        )
    elseif isempty(dist_slack) && length(ref_bus_positions) < buscount
        copyto!(PTDFm_t, BA)
        PTDFm_t[valid_ix, :] = KLU.solve!(K, PTDFm_t[valid_ix, :])
        PTDFm_t[collect(ref_bus_positions), :] .= 0.0
        return PTDFm_t
    elseif length(dist_slack) == buscount
        @info "Distributed bus"
        copyto!(PTDFm_t, BA)
        PTDFm_t[valid_ix, :] = KLU.solve!(K, PTDFm_t[valid_ix, :])
        PTDFm_t[ref_bus_positions, :] .= 0.0
        slack_array = dist_slack / sum(dist_slack)
        slack_array = reshape(slack_array, 1, buscount)
        return PTDFm_t - (PTDFm_t * slack_array) * ones(1, buscount)
    else
        error("Distributed bus specification doesn't match the number of buses.")
    end
end

"""
Computes the PTDF matrix by means of the KLU.LU factorization for sparse matrices.

# Arguments
- `branches`:
        vector of the System AC branches
- `nodes::Vector{PSY.Bus}`:
        vector of the System buses
- `bus_lookup::Dict{Int, Int}`:
        dictionary mapping the bus numbers with their enumerated indexes.
- `dist_slack::Vector{Float64}`:
        vector containing the weights for the distributed slacks.
"""
function calculate_PTDF_matrix_KLU(
    branches,
    nodes::Vector{PSY.Bus},
    bus_lookup::Dict{Int, Int},
    dist_slack::Vector{Float64})
    A, ref_bus_positions = calculate_A_matrix(branches, nodes)
    BA = calculate_BA_matrix(branches, bus_lookup)
    PTDFm = _calculate_PTDF_matrix_KLU(A, BA, ref_bus_positions, dist_slack)
    return PTDFm, A
end

"""
!!! MISSING DOCUMENTATION !!!
"""
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
Funciton for internal use only.

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
    A::Matrix{Int8},
    BA::Matrix{T},
    ref_bus_positions::Set{Int},
    dist_slack::Vector{Float64}) where {T <: Union{Float32, Float64}}
    linecount = size(BA, 2)
    buscount = size(BA, 1)
    # Use dense calculation of ABA
    valid_ixs = setdiff(1:buscount, ref_bus_positions)
    ABA = transpose(A[:, valid_ixs]) * transpose(BA[valid_ixs, :])
    # get LU factorization matrices
    (ABA, bipiv, binfo) = getrf!(ABA)
    _binfo_check(binfo)
    PTDFm_t = zeros(buscount, linecount)

    if !isempty(dist_slack) && length(ref_bus_positions) != 1
        error(
            "Distibuted slack is not supported for systems with multiple reference buses.",
        )
    elseif isempty(dist_slack) && length(ref_bus_positions) < buscount
        PTDFm_t[valid_ixs, :] = gemm(
            'N',
            'N',
            getri!(ABA, bipiv),
            BA[valid_ixs, :],
        )
        return PTDFm_t
    elseif length(dist_slack) == buscount
        @info "Distributed bus"
        slack_array = dist_slack / sum(dist_slack)
        slack_array = reshape(slack_array, buscount, 1)
        return PTDFm_t -
               gemm('N', 'N', gemm('N', 'N', PTDFm, slack_array), ones(1, buscount))
    else
        error("Distributed bus specification doesn't match the number of buses.")
    end

    return PTDFm_t
end

"""
Computes the PTDF matrix by means of the LAPACK and BLAS functions for dense matrices.

# Arguments
- `branches`:
        vector of the System AC branches
- `nodes::Vector{PSY.Bus}`:
        vector of the System buses
- `bus_lookup::Dict{Int, Int}`:
        dictionary mapping the bus numbers with their enumerated indexes.
- `dist_slack::Vector{Float64}`:
        vector containing the weights for the distributed slacks.
"""
function calculate_PTDF_matrix_DENSE(
    branches,
    nodes::Vector{PSY.Bus},
    bus_lookup::Dict{Int, Int},
    dist_slack::Vector{Float64})
    A, ref_bus_positions = calculate_A_matrix(branches, nodes)
    BA = Matrix(calculate_BA_matrix(branches, bus_lookup))
    PTDFm = _calculate_PTDF_matrix_DENSE(Matrix(A), BA, ref_bus_positions, dist_slack)
    return PTDFm, A
end

"""
Funciton for internal use only.

Computes the PTDF matrix by means of the MKL Pardiso for dense matrices.

# Arguments
- `A::SparseArrays.SparseMatrixCSC{Int8, Int}`:
        Incidence Matrix
- `BA::SparseArrays.SparseMatrixCSC{Float64, Int}`:
        BA matrix
- `ref_bus_positions::Set{Int}`:
        vector containing the indexes of the referece slack buses.
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
    # Here add the subnetwork detection
    Ix = Matrix(
        1.0I,
        buscount - length(ref_bus_positions),
        buscount - length(ref_bus_positions),
    )
    ABA_inv = zeros(Float64, size(Ix))

    ps = Pardiso.MKLPardisoSolver()
    Pardiso.set_iparm!(ps, 59, 2)
    Pardiso.set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)

    Pardiso.solve!(ps, ABA_inv, ABA, Ix)
    PTDFm_t = zeros(buscount, linecount)
    row_idx = setdiff(1:size(PTDFm_t, 1), ref_bus_positions)

    if !isempty(dist_slack) && length(ref_bus_positions) != 1
        error(
            "Distibuted slack is not supported for systems with multiple reference buses.",
        )
    elseif isempty(dist_slack) && length(ref_bus_positions) < buscount
        PTDFm_t[row_idx, :] = ABA_inv * @view BA[row_idx, :]
    elseif length(dist_slack) == buscount
        @info "Distributed bus"
        PTDFm_t[row_idx, :] = ABA_inv * @view BA[row_idx, :]
        slack_array = dist_slack / sum(dist_slack)
        slack_array = reshape(slack_array, buscount, 1)
        PTDFm_t = PTDFm_t - (PTDFm_t * slack_array) * ones(1, buscount)
    else
        error("Distributed bus specification doesn't match the number of buses.")
    end

    return PTDFm_t
end

"""
Computes the PTDF matrix by means of the MKL Pardiso for dense matrices.

# Arguments
- `branches`:
        vector of the System AC branches
- `nodes::Vector{PSY.Bus}`:
        vector of the System buses
- `bus_lookup::Dict{Int, Int}`:
        dictionary mapping the bus numbers with their enumerated indexes.
- `dist_slack::Vector{Float64}`:
        vector containing the weights for the distributed slacks.
"""
function calculate_PTDF_matrix_MKLPardiso(
    branches,
    nodes::Vector{PSY.Bus},
    bus_lookup::Dict{Int, Int},
    dist_slack::Vector{Float64})
    A, ref_bus_positions = calculate_A_matrix(branches, nodes)
    BA = calculate_BA_matrix(branches, bus_lookup)
    PTDFm = _calculate_PTDF_matrix_MKLPardiso(A, BA, ref_bus_positions, dist_slack)
    return PTDFm, A
end

"""
Builds the PTDF matrix from a group of branches and nodes. The return is a PTDF array indexed with the bus numbers.

# Arguments
- `branches`:
        vector of the System AC branches
- `nodes::Vector{PSY.Bus}`:
        vector of the System buses
- `dist_slack::Vector{Float64}`:
        vector of weights to be used as distributed slack bus.
        The distributed slack vector has to be the same length as the number of buses
- `linear_solver::String`:
        Linear solver to be used. Options are "Dense", "KLU" and "MKLPardiso
- `tol::Float64`:
        Tolerance to eliminate entries in the PTDF matrix (default eps())
"""
function PTDF(
    branches,
    nodes::Vector{PSY.Bus};
    dist_slack::Vector{Float64} = Float64[],
    linear_solver::String = "KLU",
    tol::Float64 = eps())
    validate_linear_solver(linear_solver)
    #Get axis names
    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in nodes]
    axes = (bus_ax, line_ax)
    M, bus_ax_ref = calculate_adjacency(branches, nodes)
    ref_bus_positions = find_slack_positions(nodes)
    subnetworks = find_subnetworks(M, bus_ax)
    if length(subnetworks) > 1
        @info "Network is not connected, using subnetworks"
        subnetworks = assing_reference_buses(subnetworks, ref_bus_positions)
    end
    look_up = (bus_ax_ref, make_ax_ref(line_ax))
    S, _ = _buildptdf(
        branches,
        nodes,
        look_up[1],
        dist_slack,
        linear_solver,
    )
    if tol > eps()
        return PTDF(
            sparsify(S, tol),
            axes,
            look_up,
            subnetworks,
            ref_bus_positions,
            Ref(tol),
        )
    end
    return PTDF(S, axes, look_up, subnetworks, ref_bus_positions, Ref(tol))
end

"""
Builds the PTDF matrix from a system. The return is a PTDF array indexed with the bus numbers.

# Arguments
- `sys::PSY.System`:
        Power Systems system
"""
function PTDF(
    sys::PSY.System;
    kwargs...,
)
    branches = get_ac_branches(sys)
    nodes = get_buses(sys)
    return PTDF(branches, nodes; kwargs...)
end

"""
Builds the PTDF matrix from a system. The return is a PTDF array indexed with the bus numbers.

# Arguments
- `A::IncidenceMatrix`:
        Incidence Matrix (full structure)
- `BA::BA_Matrix`:
        BA matrix (full structure)
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
    dist_slack::Vector{Float64} = Float64[],
    linear_solver = "KLU",
    tol::Float64 = eps())
    validate_linear_solver(linear_solver)
    S = _buildptdf_from_matrices(A, BA.data, dist_slack, linear_solver)
    axes = (A.axes[2], A.axes[1])
    lookup = (A.lookup[2], A.lookup[1])
    @warn "PTDF creates via other matrices doesn't compute the subnetworks"
    if tol > eps()
        return PTDF(
            sparsify(S, tol),
            axes,
            lookup,
            Dict{Int, Set{Int}}(),
            BA.ref_bus_positions,
            Ref(tol),
        )
    else
        return PTDF(S, axes, lookup, Dict{Int, Set{Int}}(), BA.ref_bus_positions, Ref(tol))
    end
end

# PTDF stores the transposed matrix. Overload indexing and how data is exported.
function Base.getindex(A::PTDF, line, bus)
    i, j = to_index(A, bus, line)
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

function get_branch_ax(ptdf::PTDF)
    return ptdf.axes[2]
end

function get_bus_ax(ptdf::PTDF)
    return ptdf.axes[1]
end
