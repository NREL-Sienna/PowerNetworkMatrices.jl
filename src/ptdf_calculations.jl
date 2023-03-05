"""
Power Transfer Distribution Factors (PTDF) indicate the incremental change in real power that occurs on transmission lines due to real power injections changes at the buses.

The PTDF struct is indexed using the Bus numbers and branch names
"""
struct PTDF{Ax, L <: NTuple{2, Dict}, T <: Real} <: PowerNetworkMatrix{Real}
    data::Array{T, 2}
    axes::Ax
    lookup::L
end

function binfo_check(binfo::Int)
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

# version 1: use same structure as already present ###########################

function _buildptdf(
    branches,
    nodes::Vector{PSY.Bus},
    dist_slack::Vector{Float64},
    linear_solver::String = "Dense")
    if linear_solver == "KLU"
        PTDFm, A = calculate_PTDF_matrix_KLU(branches, nodes, dist_slack)
    elseif linear_solver == "Dense"
        PTDFm, A = calculate_PTDF_matrix_DENSE(branches, nodes, dist_slack)
    elseif linear_solver == "MKLPardiso"
        PTDFm, A = calculate_PTDF_matrix_MKLPardiso(branches, nodes, dist_slack)
    end

    return PTDFm, A
end

function _buildptdf(
    A::IncidenceMatrix,
    BA::SparseArrays.SparseMatrixCSC{T, Int32} where {T <: Union{Float32, Float64}},
    dist_slack::Vector{Float64},
    linear_solver::String)
    if linear_solver == "KLU"
        PTDFm, A = calculate_PTDF_matrix_KLU(A, BA, dist_slack)
    elseif linear_solver == "Dense"
        PTDFm, A = calculate_PTDF_matrix_DENSE(A, BA, dist_slack)
    elseif linear_solver == "MKLPardiso"
        PTDFm, A = calculate_PTDF_matrix_MKLPardiso(A, BA, dist_slack)
    end

    return PTDFm, A
end

"""
Builds the PTDF matrix from a group of branches and nodes. The return is a PTDF array indexed with the bus numbers.

# Keyword arguments
- `dist_slack::Vector{Float64}`: Vector of weights to be used as distributed slack bus.
    The distributed slack vector has to be the same length as the number of buses
"""
function PTDF(
    branches,
    nodes::Vector{PSY.Bus},
    dist_slack::Vector{Float64} = [0.1];
    linear_solver::String = "Dense")

    #Get axis names
    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in nodes]
    S, _ = _buildptdf(branches, nodes, dist_slack, linear_solver)
    axes = (line_ax, bus_ax)
    look_up = (_make_ax_ref(line_ax), _make_ax_ref(bus_ax))

    return PTDF(S, axes, look_up)
end

"""
Builds the PTDF matrix from a system. The return is a PTDF array indexed with the bus numbers.

# Keyword arguments
- `dist_slack::Vector{Float64}`: Vector of weights to be used as distributed slack bus.
    The distributed slack vector has to be the same length as the number of buses
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
    S, _ = _buildptdf(A, BA.data, dist_slack, linear_solver)

    return PTDF(S, A.axes, A.lookup)
end
