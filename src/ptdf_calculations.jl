"""
Power Transfer Distribution Factors (PTDF) indicate the incremental change in real power that occurs on transmission lines due to real power injections changes at the buses.

The PTDF struct is indexed using the Bus numbers and branch names
"""
struct PTDF{Ax, L <: NTuple{2, Dict}, T} <: PowerNetworkMatrix{T}
    data::AbstractArray{T, 2}
    axes::Ax
    lookup::L
    tol::Float64
end

function _buildptdf(
    branches,
    nodes::Vector{PSY.Bus},
    bus_lookup::Dict{Int64, Int64},
    dist_slack::Vector{Float64},
    linear_solver::String = "Dense")
    if linear_solver == "KLU"
        PTDFm, A = calculate_PTDF_matrix_KLU(branches, nodes, bus_lookup, dist_slack)
    elseif linear_solver == "Dense"
        PTDFm, A = calculate_PTDF_matrix_DENSE(branches, nodes, bus_lookup, dist_slack)
    elseif linear_solver == "MKLPardiso"
        PTDFm, A = calculate_PTDF_matrix_MKLPardiso(branches, nodes, bus_lookup, dist_slack)
    end

    return PTDFm, A

end

function _buildptdf_from_matrices(
    A::IncidenceMatrix,
    BA::SparseArrays.SparseMatrixCSC{T, Int32} where {T <: Union{Float32, Float64}},
    dist_slack::Vector{Float64},
    linear_solver::String)
    if linear_solver == "KLU"
        PTDFm = _calculate_PTDF_matrix_KLU(A.data, BA, A.slack_positions, dist_slack)
    elseif linear_solver == "Dense"
        # Convert SparseMatrices to Dense
        PTDFm = _calculate_PTDF_matrix_DENSE(Matrix(A.data), Matrix(BA), A.slack_positions, dist_slack)
    elseif linear_solver == "MKLPardiso"
        PTDFm = _calculate_PTDF_matrix_MKLPardiso(A.data, BA, A.slack_positions, dist_slack)
    end

    return PTDFm
end

# PTDF evaluation ############################################################
function _calculate_PTDF_matrix_KLU(
    A::SparseArrays.SparseMatrixCSC{Int8, Int32},
    BA::SparseArrays.SparseMatrixCSC{Float64, Int32},
    slack_positions::Vector{Int64},
    dist_slack::Vector{Float64})

    linecount = size(BA, 1)
    buscount = size(BA, 2)

    ABA = calculate_ABA_matrix(A, BA, slack_positions)
    K = klu(ABA)
    Ix = Matrix(1.0I, size(ABA, 1), size(ABA, 1))
    ABA_inv = zeros(Float64, size(ABA, 1), size(ABA, 2))
    ldiv!(ABA_inv, K, Ix)
    PTDFm = zeros(size(BA, 1), size(BA, 2) + 1)
    slack_position = slack_positions[1]

    if dist_slack[1] == 0.1 && length(dist_slack) == 1
        PTDFm[:, setdiff(1:end, slack_positions[1])] = BA * ABA_inv
    elseif dist_slack[1] != 0.1 && length(dist_slack) == buscount
        @info "Distributed bus"
        PTDFm[:, setdiff(1:end, slack_positions[1])] = BA * ABA_inv
        slack_array = dist_slack / sum(dist_slack)
        slack_array = reshape(slack_array, buscount, 1)
        PTDFm = PTDFm - (PTDFm * slack_array) * ones(1, buscount)
    elseif length(slack_position) == 0
        @warn("Slack bus not identified in the Bus/Nodes list, can't build PTDF")
        PTDFm = Array{Float64, 2}(undef, linecount, buscount)
    else
        @assert false
    end

    return PTDFm
end

function calculate_PTDF_matrix_KLU(
    branches,
    nodes::Vector{PSY.Bus},
    bus_lookup::Dict{Int64, Int64},
    dist_slack::Vector{Float64})
    A, slack_positions = calculate_A_matrix(branches, nodes)
    BA = calculate_BA_matrix(branches, slack_positions, bus_lookup)
    PTDFm = _calculate_PTDF_matrix_KLU(A, BA, slack_positions, dist_slack)
    return PTDFm, A
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

function _calculate_PTDF_matrix_DENSE(
    A::Matrix{Int8},
    BA::Matrix{T},
    slack_position::Vector{Int64},
    dist_slack::Vector{Float64}) where {T <: Union{Float32, Float64}}

    # Use dense calculation of ABA
    ABA = A[:, setdiff(1:end, slack_position[1])]' * BA
    linecount = size(BA, 1)
    buscount = size(BA, 2)
    # get LU factorization matrices
    if dist_slack[1] == 0.1 && length(dist_slack) == 1
        (ABA, bipiv, binfo) = getrf!(ABA)
        _binfo_check(binfo)
        PTDFm = gemm(
            'N',
            'N',
            BA,
            getri!(ABA, bipiv),
        )
        @views PTDFm =
            hcat(
                PTDFm[:, 1:(slack_position[1] - 1)],
                zeros(linecount),
                PTDFm[:, slack_position[1]:end],
            )
    elseif dist_slack[1] != 0.1 && length(dist_slack) == buscount
        @info "Distributed bus"
        (ABA, bipiv, binfo) = getrf!(ABA)
        _binfo_check(binfo)
        PTDFm = gemm(
            'N',
            'N',
            BA,
            getri!(ABA, bipiv),
        )
        @views PTDFm =
            hcat(
                PTDFm[:, 1:(slack_position[1] - 1)],
                zeros(linecount),
                PTDFm[:, slack_position[1]:end],
            )
        slack_array = dist_slack / sum(dist_slack)
        slack_array = reshape(slack_array, buscount, 1)
        PTDFm =
            PTDFm - gemm('N', 'N', gemm('N', 'N', PTDFm, slack_array), ones(1, buscount))
    elseif length(slack_position) == 0
        @warn("Slack bus not identified in the Bus/Nodes list, can't build PTDF")
        PTDFm = Array{Float64, 2}(undef, linecount, buscount)
    else
        @assert false
    end

    return PTDFm
end

function calculate_PTDF_matrix_DENSE(
    branches,
    nodes::Vector{PSY.Bus},
    bus_lookup::Dict{Int64, Int64},
    dist_slack::Vector{Float64})
    A, slack_positions = calculate_A_matrix(branches, nodes)
    BA = Matrix(calculate_BA_matrix(branches, slack_positions, bus_lookup))
    PTDFm = _calculate_PTDF_matrix_DENSE(Matrix(A), BA, slack_positions, dist_slack)
    return PTDFm, A
end

function _calculate_PTDF_matrix_MKLPardiso(
    A::SparseArrays.SparseMatrixCSC{Int8, Int32},
    BA::SparseArrays.SparseMatrixCSC{Float64, Int32},
    slack_positions::Vector{Int64},
    dist_slack::Vector{Float64})

    ps = Pardiso.MKLPardisoSolver()

    linecount = size(BA, 1)
    buscount = size(BA, 2)

    ABA = calculate_ABA_matrix(A, BA, slack_positions)
    Ix = Matrix(1.0I, size(ABA, 1), size(ABA, 1))
    ABA_inv = zeros(Float64, size(ABA, 1), size(ABA, 2))
    Pardiso.solve!(ps, ABA_inv, ABA, Ix)
    PTDFm = zeros(size(BA, 1), size(BA, 2) + 1)
    slack_position = slack_positions[1]

    if dist_slack[1] == 0.1 && length(dist_slack) == 1
        PTDFm[:, setdiff(1:end, slack_positions[1])] = BA * ABA_inv
    elseif dist_slack[1] != 0.1 && length(dist_slack) == buscount
        @info "Distributed bus"
        PTDFm[:, setdiff(1:end, slack_positions[1])] = BA * ABA_inv
        slack_array = dist_slack / sum(dist_slack)
        slack_array = reshape(slack_array, buscount, 1)
        PTDFm = PTDFm - (PTDFm * slack_array) * ones(1, buscount)
    elseif length(slack_position) == 0
        @warn("Slack bus not identified in the Bus/Nodes list, can't build PTDF")
        PTDFm = Array{Float64, 2}(undef, linecount, buscount)
    else
        @assert false
    end

    return PTDFm
end

function calculate_PTDF_matrix_MKLPardiso(
    branches,
    nodes::Vector{PSY.Bus},
    bus_lookup::Dict{Int64, Int64},
    dist_slack::Vector{Float64})

    A, slack_positions = calculate_A_matrix(branches, nodes)
    BA = calculate_BA_matrix(branches, slack_positions, bus_lookup)
    PTDFm = _calculate_PTDF_matrix_MKLPardiso(A, BA, slack_positions, dist_slack)
    return PTDFm, A
end

"""
Builds the PTDF matrix from a group of branches and nodes. The return is a PTDF array indexed with the bus numbers.

# Keyword arguments
- `dist_slack::Vector{Float64}`: Vector of weights to be used as distributed slack bus.
    The distributed slack vector has to be the same length as the number of buses
- `linear_solver::String`: Linear solver to be used. Options are "Dense", "KLU" and "MKLPardiso
- `tol::Float64`: Tolerance to eliminate entries in the PTDF matrix (default eps())
"""
function PTDF(
    branches,
    nodes::Vector{PSY.Bus},
    dist_slack::Vector{Float64} = [0.1];
    linear_solver::String = "Dense")

    #Get axis names
    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in nodes]
    axes = (line_ax, bus_ax)
    look_up = (_make_ax_ref(line_ax), _make_ax_ref(bus_ax))
    S, _ = _buildptdf(branches, nodes, look_up[2], dist_slack, linear_solver)
    return PTDF(S, axes, look_up, eps())
end

"""
Builds the PTDF matrix from a system. The return is a PTDF array indexed with the bus numbers.

# Keyword arguments
- `dist_slack::Vector{Float64}`: Vector of weights to be used as distributed slack bus.
    The distributed slack vector has to be the same length as the number of buses
- `linear_solver::String`: Linear solver to be used. Options are "Dense", "KLU" and "MKLPardiso
- `tol::Float64`: Tolerance to eliminate entries in the PTDF matrix (default eps())
"""
function PTDF(
    sys::PSY.System,
    dist_slack::Vector{Float64} = [0.1];
    linear_solver::String = "Dense")
    branches = get_ac_branches(sys)
    nodes = get_buses(sys)
    validate_linear_solver(linear_solver)

    return PTDF(branches, nodes, dist_slack; linear_solver)
end

# version 2: use BA and ABA fucntions created before #########################

function PTDF(
    A::IncidenceMatrix,
    BA::BA_Matrix,
    dist_slack::Vector{Float64} = [0.1];
    linear_solver = "Dense")
    validate_linear_solver(linear_solver)
    S = _buildptdf_from_matrices(A, BA.data, dist_slack, linear_solver)

    return PTDF(S, A.axes, A.lookup, eps())
end
