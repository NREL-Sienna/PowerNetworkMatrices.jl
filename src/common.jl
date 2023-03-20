# get branches from system
function get_ac_branches(sys::PSY.System)
    # Filter out DC Branches here
    return sort!(
        collect(PSY.get_components(PSY.get_available, PSY.ACBranch, sys));
        by = x -> (PSY.get_number(PSY.get_arc(x).from), PSY.get_number(PSY.get_arc(x).to)),
    )
end

# get buses from system
function get_buses(sys::PSY.System)::Vector{PSY.Bus}
    return sort!(collect(PSY.get_components(PSY.Bus, sys)); by = x -> PSY.get_number(x))
end

# get slack bus
function find_slack_positions(nodes)
    bus_lookup = make_ax_ref(nodes)
    slack_position = sort([
        bus_lookup[PSY.get_number(n)]
        for n in nodes if PSY.get_bustype(n) == BusTypes.REF
    ])
    if length(slack_position) == 0
        error("Slack bus not identified in the Bus/Nodes list, can't build NetworkMatrix")
    end
    return slack_position
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
    ref_bus_positions = find_slack_positions(nodes)

    A_I = Int[]
    A_J = Int[]
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

    return SparseArrays.sparse(A_I, A_J, A_V), ref_bus_positions
end

function calculate_adjacency(branches, nodes::Vector{PSY.Bus})
    bus_ax = PSY.get_number.(nodes)
    return calculate_adjacency(branches, nodes, make_ax_ref(bus_ax))
end

function calculate_adjacency(
    branches,
    nodes::Vector{PSY.Bus},
    bus_lookup::Dict{Int, Int},
)
    buscount = length(nodes)
    a = SparseArrays.spzeros(Int8, buscount, buscount)

    for b in branches
        (fr_b, to_b) = get_bus_indices(b, bus_lookup)
        a[fr_b, to_b] = 1
        a[to_b, fr_b] = -1
        a[fr_b, fr_b] = 1
        a[to_b, to_b] = 1
    end

    # Return both for type stability
    return a, bus_lookup
end

# BA matrix evaluation #######################################################
function calculate_BA_matrix(
    branches,
    ref_bus_positions::Vector{Int},
    bus_lookup::Dict{Int, Int})
    BA_I = Int[]
    BA_J = Int[]
    BA_V = Float64[]

    for (ix, b) in enumerate(branches)
        if isa(b, PSY.DCBranch)
            @warn("PTDF construction ignores DC-Lines")
            continue
        end

        (fr_b, to_b) = get_bus_indices(b, bus_lookup)
        b_val = PSY.get_series_susceptance(b)

        if fr_b ∉ ref_bus_positions
            check_ = sum(fr_b .> ref_bus_positions)
            push!(BA_I, ix)
            push!(BA_J, fr_b - check_)
            push!(BA_V, b_val)
        end

        if to_b ∉ ref_bus_positions
            check_ = sum(to_b .> ref_bus_positions)
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
    A::SparseArrays.SparseMatrixCSC{Int8, Int},
    BA::SparseArrays.SparseMatrixCSC{T, Int} where {T <: Union{Float32, Float64}},
    ref_bus_positions::Vector{Int})
    return A[:, setdiff(1:end, ref_bus_positions)]' * BA
end
