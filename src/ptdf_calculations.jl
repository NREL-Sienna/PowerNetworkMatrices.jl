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
    ref_bus_positions::Vector{Int}
    tol::Base.RefValue{Float64}
end

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
Makes the PTDF matrix sparse given a certain tolerance "tol".
"""
function make_sparse_PTDF(mat::PTDF{Ax, L, Matrix{Float64}}, tol::Float64) where {Ax, L}
    new_mat = sparsify(mat.data, tol)
    return PTDF(new_mat, mat.axes, mat.lookup, mat.subnetworks, Ref(tol))
end

function _buildptdf(
    branches,
    nodes::Vector{PSY.Bus},
    bus_lookup::Dict{Int, Int},
    dist_slack::Vector{Float64},
    ref_bus_positions::Vector{Int},
    linear_solver::String)
    if linear_solver == "KLU"
        PTDFm, A = calculate_PTDF_matrix_KLU(
            branches,
            nodes,
            bus_lookup,
            ref_bus_positions,
            dist_slack,
        )
    elseif linear_solver == "Dense"
        PTDFm, A = calculate_PTDF_matrix_DENSE(
            branches,
            nodes,
            bus_lookup,
            ref_bus_positions,
            dist_slack,
        )
    elseif linear_solver == "MKLPardiso"
        PTDFm, A = calculate_PTDF_matrix_MKLPardiso(
            branches,
            nodes,
            bus_lookup,
            ref_bus_positions,
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
- `ref_bus_positions::Vector{Int}`:
        vector containing the indexes of the reference slack buses.
- `dist_slack::Vector{Float64}`:
        vector containing the weights for the distributed slacks.
"""
function _calculate_PTDF_matrix_KLU(
    A::SparseArrays.SparseMatrixCSC{Int8, Int},
    BA::SparseArrays.SparseMatrixCSC{Float64, Int},
    ref_bus_positions::Vector{Int},
    dist_slack::Vector{Float64})
    linecount = size(BA, 1)
    buscount = size(BA, 2)

    ABA = calculate_ABA_matrix(A, BA, ref_bus_positions)
    # Here add the subnetwork detection
    K = klu(ABA)
    Ix = Matrix(
        1.0I,
        buscount - length(ref_bus_positions),
        buscount - length(ref_bus_positions),
    )
    ABA_inv = zeros(Float64, size(Ix))
    ldiv!(ABA_inv, K, Ix)
    PTDFm = zeros(linecount, buscount)
    if !isempty(dist_slack) && length(ref_bus_positions) != 1
        error(
            "Distibuted slack is not supported for systems with multiple reference buses.",
        )
    elseif isempty(dist_slack) && length(ref_bus_positions) < buscount
        A_mul_B!(PTDFm, BA, ABA_inv, ref_bus_positions)
    elseif length(dist_slack) == buscount
        @info "Distributed bus"
        A_mul_B!(PTDFm, BA, ABA_inv, ref_bus_positions)
        slack_array = dist_slack / sum(dist_slack)
        slack_array = reshape(slack_array, buscount, 1)
        PTDFm = PTDFm - (PTDFm * slack_array) * ones(1, buscount)
    else
        error("Distributed bus specification doesn't match the number of buses.")
    end
    return PTDFm
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
    ref_bus_positions::Vector{Int},
    dist_slack::Vector{Float64})
    A, ref_bus_positions = calculate_A_matrix(branches, nodes, ref_bus_positions)
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
- `ref_bus_positions::Vector{Int}`:
        vector containing the indexes of the reference slack buses.
- `dist_slack::Vector{Float64})`:
        vector containing the weights for the distributed slacks.
"""
function _calculate_PTDF_matrix_DENSE(
    A::Matrix{Int8},
    BA::Matrix{T},
    ref_bus_positions::Vector{Int},
    dist_slack::Vector{Float64}) where {T <: Union{Float32, Float64}}

    linecount = size(BA, 1)
    buscount = size(BA, 2)
    # Use dense calculation of ABA
    valid_ixs = setdiff(1:buscount, ref_bus_positions)
    ABA = zeros(Float64, length(valid_ixs), length(valid_ixs))
    ABA[:] = (transpose(A) * BA)[valid_ixs, valid_ixs]
    # get LU factorization matrices
    if !isempty(dist_slack) && length(ref_bus_positions) != 1
        error(
            "Distibuted slack is not supported for systems with multiple reference buses.",
        )
    end

    (ABA, bipiv, binfo) = getrf!(ABA)
    _binfo_check(binfo)
    PTDFm = zeros(Float64, linecount, buscount)
    PTDFm[:, valid_ixs] = gemm(
        'N',
        'N',
        BA[:, valid_ixs],
        getri!(ABA, bipiv),
    )

    if length(dist_slack) == buscount
        @info "Distributed bus"
        slack_array = dist_slack / sum(dist_slack)
        slack_array = reshape(slack_array, buscount, 1)
        PTDFm =
            PTDFm - gemm('N', 'N', gemm('N', 'N', PTDFm, slack_array), ones(1, buscount))
    end

    return PTDFm
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
    ref_bus_positions::Vector{Int},
    dist_slack::Vector{Float64})
    A, ref_bus_positions = calculate_A_matrix(branches, nodes, ref_bus_positions)
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
- `ref_bus_positions::Vector{Int}`:
        vector containing the indexes of the referece slack buses.
- `dist_slack::Vector{Float64}`:
        vector containing the weights for the distributed slacks.
"""
function _calculate_PTDF_matrix_MKLPardiso(
    A::SparseArrays.SparseMatrixCSC{Int8, Int},
    BA::SparseArrays.SparseMatrixCSC{Float64, Int},
    ref_bus_positions::Vector{Int},
    dist_slack::Vector{Float64})
    ps = Pardiso.MKLPardisoSolver()

    linecount = size(BA, 1)
    buscount = size(BA, 2)

    ABA = calculate_ABA_matrix(A, BA, ref_bus_positions)
    # Here add the subnetwork detection
    Ix = Matrix(
        1.0I,
        buscount - length(ref_bus_positions),
        buscount - length(ref_bus_positions),
    )
    ABA_inv = zeros(Float64, size(Ix))
    Pardiso.solve!(ps, ABA_inv, ABA, Ix)
    PTDFm = zeros(linecount, buscount)

    if !isempty(dist_slack) && length(ref_bus_positions) != 1
        error(
            "Distibuted slack is not supported for systems with multiple reference buses.",
        )
    elseif isempty(dist_slack) && length(ref_bus_positions) < buscount
        A_mul_B!(PTDFm, BA, ABA_inv, ref_bus_positions)
    elseif length(dist_slack) == buscount
        @info "Distributed bus"
        A_mul_B!(PTDFm, BA, ABA_inv, ref_bus_positions)
        slack_array = dist_slack / sum(dist_slack)
        slack_array = reshape(slack_array, buscount, 1)
        PTDFm = PTDFm - (PTDFm * slack_array) * ones(1, buscount)
    else
        error("Distributed bus specification doesn't match the number of buses.")
    end

    return PTDFm
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
    ref_bus_positions::Vector{Int},
    dist_slack::Vector{Float64})
    A, ref_bus_positions = calculate_A_matrix(branches, nodes, ref_bus_positions)
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
    axes = (line_ax, bus_ax)
    M, bus_ax_ref = calculate_adjacency(branches, nodes)
    ref_bus_positions = find_slack_positions(nodes)
    subnetworks = find_subnetworks(M, bus_ax)
    if length(subnetworks) > 1
        @info "Network is not connected, using subnetworks"
        subnetworks = assing_reference_buses(subnetworks, ref_bus_positions)
    end
    look_up = (make_ax_ref(line_ax), bus_ax_ref)
    S, _ = _buildptdf(
        branches,
        nodes,
        look_up[2],
        dist_slack,
        ref_bus_positions,
        linear_solver,
    )
    if tol > eps()
        return PTDF(sparsify(S, tol), axes, look_up, Ref(tol))
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
    if tol > eps()
        return PTDF(sparsify(S, tol), axes, look_up, tol)
    end
    @warn "PTDF creates via other matrices doesn't compute the subnetworks"
    return PTDF(S, A.axes, A.lookup, Dict{Int, Set{Int}}(), BA.ref_bus_positions, Ref(tol))
end
