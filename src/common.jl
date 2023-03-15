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
    bus_lookup = make_ax_ref(nodes)
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
    bus_lookup = make_ax_ref(nodes)
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
