# get branches from system
function get_ac_branches(sys::PSY.System)
    # Filter out DC Branches here
    return sort!(
        collect(PSY.get_components(PSY.ACBranch, sys));
        by = x -> (PSY.get_number(PSY.get_arc(x).from), PSY.get_number(PSY.get_arc(x).to)),
    )
end

# get buses from system
function get_buses(sys::PSY.System)
    return sort!(collect(PSY.get_components(PSY.Bus, sys)); by = x -> PSY.get_number(x))
end

# get slack bus
function find_slack_positions(nodes)
    bus_lookup = _make_ax_ref(nodes)
    return sort([
        bus_lookup[PSY.get_number(n)]
        for n in nodes if PSY.get_bustype(n) == BusTypes.REF
    ])
end

# validate solver so be used for
function validate_linear_solver(linear_solver::String)
    if linear_solver ∉ SUPPORTED_LINEAR_SOLVERS
        error(
            "Invalid linear solver. Supported linear solvers are: $(SUPPORTED_LINEAR_SOLVERS)",
        )
    end
    return
end

# incidence matrix (A) evaluation ############################################
function calculate_A_matrix(branches, nodes::Vector{PSY.Bus})
    bus_lookup = _make_ax_ref(nodes)
    slack_positions = find_slack_positions(nodes)

    A_I = Int32[]
    A_J = Int32[]
    A_V = Int8[]

    # build incidence matrix A (lines x buses)
    for (ix, b) in enumerate(branches)
        (fr_b, to_b) = get_bus_indices(b, bus_lookup)

        # change column number
        push!(A_I, ix)
        push!(A_J, fr_b)
        push!(A_V, 1)

        push!(A_I, ix)
        push!(A_J, to_b)
        push!(A_V, -1)
    end

    return SparseArrays.sparse(A_I, A_J, A_V), slack_positions
end

# BA matrix evaluation #######################################################
function calculate_BA_matrix(
    branches,
    slack_positions::Vector{Int64},
    bus_lookup::Dict{Int64, Int64})
    BA_I = Int32[]
    BA_J = Int32[]
    BA_V = Float64[]

    for (ix, b) in enumerate(branches)
        if isa(b, PSY.DCBranch)
            @warn("PTDF construction ignores DC-Lines")
            continue
        end

        (fr_b, to_b) = get_bus_indices(b, bus_lookup)
        b_val = PSY.get_series_susceptance(b)

        if fr_b ∉ slack_positions[1]
            check_ = sum(fr_b .> slack_positions[1])
            push!(BA_I, ix)
            push!(BA_J, fr_b - check_)
            push!(BA_V, b_val)
        end

        if to_b ∉ slack_positions[1]
            check_ = sum(to_b .> slack_positions[1])
            push!(BA_I, ix)
            push!(BA_J, to_b - check_)
            push!(BA_V, -b_val)
        end
    end

    BA = SparseArrays.sparse(BA_I, BA_J, BA_V)

    return BA
end

# ABA matrix evaluation ######################################################
function calculate_ABA_matrix(
    A::SparseArrays.SparseMatrixCSC{Int8, Int32},
    BA::SparseArrays.SparseMatrixCSC{T, Int32} where {T <: Union{Float32, Float64}},
    slack_positions::Vector{Int64})
    return A[:, setdiff(1:end, slack_positions[1])]' * BA
end

# PTDF evaluation ############################################################
function calculate_PTDF_matrix_KLU(
    branches,
    nodes::Vector{PSY.Bus},
    dist_slack::Vector{Float64})
    buscount = length(nodes)
    linecount = length(branches)
    A, slack_positions = calculate_A_matrix(branches, nodes)
    bus_lookup = _make_ax_ref(nodes)

    BA = calculate_BA_matrix(branches, slack_positions, bus_lookup)
    ABA = calculate_ABA_matrix(A, BA, slack_positions)
    K = klu(ABA)
    Ix = Matrix(1.0LinearAlgebra.I, size(ABA, 1), size(ABA, 1))
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

    return PTDFm, A
end

function calculate_PTDF_matrix_KLU(
    A::IncidenceMatrix,
    BA::SparseArrays.SparseMatrixCSC{T, Int32} where {T <: Union{Float32, Float64}},
    dist_slack::Vector{Float64})
    linecount = length(A.axes[1])
    buscount = length(A.axes[2])
    slack_positions = A.slack_positions

    ABA = calculate_ABA_matrix(A.data, BA, slack_positions)
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

    return PTDFm, A.data
end

function calculate_PTDF_matrix_DENSE(
    branches,
    nodes::Vector{PSY.Bus},
    dist_slack::Vector{Float64})
    buscount = length(nodes)
    linecount = length(branches)
    bus_lookup = _make_ax_ref(nodes)
    A = zeros(Float64, buscount, linecount)
    inv_X = zeros(Float64, linecount, linecount)

    #build incidence matrix
    for (ix, b) in enumerate(branches)
        if isa(b, PSY.DCBranch)
            @warn("PTDF construction ignores DC-Lines")
            continue
        end

        (fr_b, to_b) = get_bus_indices(b, bus_lookup)
        A[fr_b, ix] = 1
        A[to_b, ix] = -1

        inv_X[ix, ix] = PSY.get_series_susceptance(b)
    end

    slacks =
        [bus_lookup[PSY.get_number(n)] for n in nodes if PSY.get_bustype(n) == BusTypes.REF]
    isempty(slacks) && error("System must have a reference bus!")

    slack_position = slacks[1]
    B = gemm(
        'N',
        'T',
        gemm('N', 'N', A[setdiff(1:end, slack_position), 1:end], inv_X),
        A[setdiff(1:end, slack_position), 1:end],
    )
    # get LU factorization matrices
    if dist_slack[1] == 0.1 && length(dist_slack) == 1
        (B, bipiv, binfo) = getrf!(B)
        binfo_check(binfo)
        PTDFm = gemm(
            'N',
            'N',
            gemm('N', 'T', inv_X, A[setdiff(1:end, slack_position), :]),
            getri!(B, bipiv),
        )
        @views PTDFm =
            hcat(
                PTDFm[:, 1:(slack_position - 1)],
                zeros(linecount),
                PTDFm[:, slack_position:end],
            )
    elseif dist_slack[1] != 0.1 && length(dist_slack) == buscount
        @info "Distributed bus"
        (B, bipiv, binfo) = getrf!(B)
        binfo_check(binfo)
        PTDFm = gemm(
            'N',
            'N',
            gemm('N', 'T', inv_X, A[setdiff(1:end, slack_position), :]),
            getri!(B, bipiv),
        )
        @views PTDFm =
            hcat(
                PTDFm[:, 1:(slack_position - 1)],
                zeros(linecount),
                PTDFm[:, slack_position:end],
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

    return PTDFm, A
end

function calculate_PTDF_matrix_MKLPardiso(
    branches,
    nodes::Vector{PSY.Bus},
    dist_slack::Vector{Float64})

    # initialize solver and number of threads to use
    ENV["OPENBLAS_NUM_THREADS"] = 10
    ENV["MKL_NUM_THREADS"] = 10
    ps = Pardiso.MKLPardisoSolver()
    set_num_threads(10)
    Pardiso.set_nprocs!(ps, 10)

    buscount = length(nodes)
    linecount = length(branches)
    A, slack_positions = calculate_A_matrix(branches, nodes)
    bus_lookup = _make_ax_ref(nodes)

    BA = calculate_BA_matrix(branches, slack_positions, bus_lookup)
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

    return PTDFm, A
end

function calculate_PTDF_matrix_MKLPardiso(
    A::IncidenceMatrix,
    BA::SparseArrays.SparseMatrixCSC{T, Int32} where {T <: Union{Float32, Float64}},
    dist_slack::Vector{Float64})

    # initialize solver and number of threads to use
    ENV["OPENBLAS_NUM_THREADS"] = 10
    ENV["MKL_NUM_THREADS"] = 10
    ps = Pardiso.MKLPardisoSolver()
    set_num_threads(10)
    Pardiso.set_nprocs!(ps, 10)

    linecount = length(A.axes[1])
    buscount = length(A.axes[2])
    slack_positions = A.slack_positions

    ABA = calculate_ABA_matrix(A.data, BA, slack_positions)
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

    return PTDFm, A.data
end
