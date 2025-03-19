function _add_to_collection!(collection::Vector{PSY.ACBranch}, branch::PSY.ACBranch)
    push!(collection, branch)
    return
end

function _add_to_collection!(
    collection::Vector{PSY.Transformer3W},
    transformer_3W::PSY.Transformer3W,
)
    push!(collection, transformer_3W)
    return
end

function _add_to_collection!(
    ::Vector{PSY.ACBranch},
    ::Union{PSY.TwoTerminalGenericHVDCLine, PSY.TwoTerminalVSCLine, PSY.TwoTerminalLCCLine},
)
    return
end

"""
Gets the AC branches from a given Systems.
"""
function get_ac_branches(
    sys::PSY.System,
    radial_branches::Set{String} = Set{String}(),
)::Vector{PSY.ACBranch}
    collection = Vector{PSY.ACBranch}()
    for br in PSY.get_components(
        x -> PSY.get_available(x) && !(typeof(x) <: PSY.Transformer3W),
        PSY.ACBranch,
        sys,
    )
        arc = PSY.get_arc(br)
        if PSY.get_bustype(arc.from) == ACBusTypes.ISOLATED
            throw(
                IS.ConflictingInputsError(
                    "Branch $(PSY.get_name(br)) is set available and connected to isolated bus $(PSY.get_name(arc.from))",
                ),
            )
        end
        if PSY.get_bustype(arc.to) == ACBusTypes.ISOLATED
            throw(
                IS.ConflictingInputsError(
                    "Branch $(PSY.get_name(br)) is set available and connected to isolated bus $(PSY.get_name(arc.to))",
                ),
            )
        end
        if PSY.get_name(br) ∉ radial_branches
            _add_to_collection!(collection, br)
        end
    end
    return sort!(collection;
        by = x -> (PSY.get_number(PSY.get_arc(x).from), PSY.get_number(PSY.get_arc(x).to)),
    )
end

"""
Gets the AC branches from a given Systems.
"""
function get_transformers_3w(
    sys::PSY.System,
)::Vector{PSY.Transformer3W}
    collection = Vector{PSY.Transformer3W}()
    for br in PSY.get_components(x -> PSY.get_available(x), PSY.Transformer3W, sys)
        ps_arc = PSY.get_primary_secondary_arc(br)
        st_arc = PSY.get_secondary_tertiary_arc(br)
        if PSY.get_bustype(ps_arc.from) == ACBusTypes.ISOLATED
            throw(
                IS.ConflictingInputsError(
                    "Branch $(PSY.get_name(br)) is set available and connected to isolated bus $(PSY.get_name(ps_arc.from))",
                ),
            )
        end
        if PSY.get_bustype(ps_arc.to) == ACBusTypes.ISOLATED
            throw(
                IS.ConflictingInputsError(
                    "Branch $(PSY.get_name(br)) is set available and connected to isolated bus $(PSY.get_name(ps_arc.to))",
                ),
            )
        end
        if PSY.get_bustype(st_arc.to) == ACBusTypes.ISOLATED
            throw(
                IS.ConflictingInputsError(
                    "Branch $(PSY.get_name(br)) is set available and connected to isolated bus $(PSY.get_name(st_arc.to))",
                ),
            )
        end
        _add_to_collection!(collection, br)
    end
    return sort!(collection;
        by = x -> (
            PSY.get_number(PSY.get_primary_secondary_arc(x).from),
            PSY.get_number(PSY.get_primary_secondary_arc(x).to),
            PSY.get_number(PSY.get_primary_tertiary_arc(x).to),
        ),
    )
end

function _next_branch_number(::ACBranch, branch_number::Int)
    return branch_number + 1
end

function _next_branch_number(::Transformer3W, branch_number::Int)
    return branch_number + 3
end

function _add_branch_to_lookup!(
    branch_lookup::Dict{String, Int},
    ::Dict{String, Vector{Int}},
    branch_type::Vector{DataType},
    branch::PSY.ACBranch,
    branch_number::Int,
)
    branch_lookup[PSY.get_name(branch)] = branch_number
    push!(branch_type, typeof(branch))
    return
end

function _add_branch_to_lookup!(
    branch_lookup::Dict{String, Int},
    transformer_3w_lookup::Dict{String, Vector{Int}},
    branch_type::Vector{DataType},
    branch::PSY.Transformer3W,
    branch_number::Int,
)
    tr3w_name = PSY.get_name(branch)
    transformer_3w_lookup[tr3w_name] = Vector{Int}(undef, 3)
    for (i, side) in enumerate(["primary", "secondary", "tertiary"])
        side_name = "$(tr3w_name)__$side"
        branch_lookup[side_name] = branch_number - 3 + i
        transformer_3w_lookup[tr3w_name] = side_name
        push!(branch_type, typeof(branch))
    end
    return
end

function get_branch_lookups(branches)
    branch_lookup = Dict{String, Int}()
    transformer_3w_lookup = Dict{String, Vector{Int}}()
    branch_type = Vector{DataType}()
    branch_number = 0
    for b in branches
        branch_number = _next_branch_number(b, branch_number)
        _add_branch_to_lookup!(
            branch_lookup,
            transformer_3w_lookup,
            branch_type,
            b,
            branch_number,
        )
    end
    return branch_lookup, transformer_3w_lookup, branch_type
end

"""
Gets the non-isolated buses from a given System
"""
function get_buses(
    sys::PSY.System,
    bus_reduction_map::Dict{Int64, Set{Int64}} = Dict{Int64, Set{Int64}}(),
)::Vector{PSY.ACBus}
    leaf_buses = Set{PSY.Int64}()
    if !isempty(bus_reduction_map)
        for vals in values(bus_reduction_map)
            union!(leaf_buses, vals)
        end
    end

    count_i = 1
    all_buses = PSY.get_components(PSY.ACBus, sys)
    buses = Vector{PSY.ACBus}(undef, length(all_buses))
    for b in all_buses
        if PSY.get_bustype(b) == ACBusTypes.ISOLATED
            continue
        end

        if PSY.get_number(b) ∈ leaf_buses
            continue
        end
        buses[count_i] = b
        count_i += 1
    end

    return sort!(deleteat!(buses, count_i:length(buses)); by = x -> PSY.get_number(x))
end

"""
Gets the indices  of the reference (slack) buses.
NOTE:
- the indices  corresponds to the columns of zeros belonging to the PTDF matrix.
- BA and ABA matrix miss the columns related to the reference buses.
"""
function find_slack_positions(buses)
    return find_slack_positions(buses, make_ax_ref(buses))
end

function find_slack_positions(buses, bus_lookup::Dict{Int, Int})::Set{Int}
    slack_position = sort([
        bus_lookup[PSY.get_number(n)]
        for n in buses if PSY.get_bustype(n) == ACBusTypes.REF
    ])
    if length(slack_position) == 0
        error("Slack bus not identified in the Bus/buses list, can't build NetworkMatrix")
    end
    return Set{Int}(slack_position)
end

"""
Validates if the selected linear solver is supported.
"""
function validate_linear_solver(linear_solver::String)
    if linear_solver ∉ SUPPORTED_LINEAR_SOLVERS
        error(
            "Invalid linear solver. Supported linear solvers are: $(SUPPORTED_LINEAR_SOLVERS)",
        )
    end
    return
end

"""
Evaluates the Incidence matrix A given the branches and node of a System.

# Arguments
- `branches`:
        vector containing the branches of the considered system (should be AC branches).
- `buses::Vector{PSY.ACBus}`:
        vector containing the buses of the considered system.

NOTE:
- the matrix features all the columns, including the ones related to the
  reference buses (each column is related to a system's bus).
"""
function calculate_A_matrix(
    branches,
    buses::Vector{PSY.ACBus},
)
    ref_bus_positions = find_slack_positions(buses)
    bus_lookup = make_ax_ref(buses)

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

"""
Evaluates the Adjacency matrix given the banches and buses of a given System.

# Arguments
- `branches`:
        vector containing the branches of the considered system (should be AC branches).
- `buses::Vector{PSY.ACBus}`:
        vector containing the buses of the considered system.
"""
function calculate_adjacency(branches, buses::Vector{PSY.ACBus})
    bus_ax = PSY.get_number.(buses)
    return calculate_adjacency(branches, buses, make_ax_ref(bus_ax))
end

"""
Evaluates the Adjacency matrix given the System's banches, buses and bus_lookup.

NOTE:
- bus_lookup is a dictionary mapping the bus numbers (as shown in the Systems)
  with their enumerated indxes.
"""
function calculate_adjacency(
    branches,
    buses::Vector{PSY.ACBus},
    bus_lookup::Dict{Int, Int},
)
    buscount = length(buses)
    a = SparseArrays.spzeros(Int8, buscount, buscount)

    for b in branches
        fr_b, to_b = get_bus_indices(b, bus_lookup)
        a[fr_b, to_b] = 1
        a[to_b, fr_b] = -1
    end

    # If a line is disconnected needs to check for the buses correctly
    for i in 1:buscount
        if PSY.get_bustype(buses[i]) == ACBusTypes.ISOLATED
            continue
        end
        a[i, i] = 1
    end

    # Return both for type stability
    return a, bus_lookup
end

"""
Evaluates the transposed BA matrix given the System's banches, reference bus positions and bus_lookup.

# Arguments
- `branches`:
        vector containing the branches of the considered system (should be AC branches).
- `ref_bus_positions::Set{Int}`:
        Vector containing the indexes of the columns of the BA matrix corresponding
        to the reference buses
- `bus_lookup::Dict{Int, Int}`:
        dictionary mapping the bus numbers with their enumerated indexes.
"""
function calculate_BA_matrix(
    branches,
    bus_lookup::Dict{Int, Int})
    BA_I = Int[]
    BA_J = Int[]
    BA_V = Float64[]

    for (ix, b) in enumerate(branches)
        (fr_b, to_b) = get_bus_indices(b, bus_lookup)
        b_val = PSY.get_series_susceptance(b)

        if !isfinite(b_val)
            error("Invalid value for branch $(PSY.summary(b)), $b_val")
        end

        push!(BA_I, fr_b)
        push!(BA_J, ix)
        push!(BA_V, b_val)

        push!(BA_I, to_b)
        push!(BA_J, ix)
        push!(BA_V, -b_val)
    end

    BA = SparseArrays.sparse(BA_I, BA_J, BA_V)

    return BA
end

"""
Evaluates the ABA matrix given the System's Incidence matrix (A), BA matrix and
reference bus positions.

# Arguments
- `A::SparseArrays.SparseMatrixCSC{Int8, Int}`:
        Incidence matrix.
- `BA::SparseArrays.SparseMatrixCSC{Float64, Int}`
        BA matrix.

NOTE:
- evaluates A with "calculate_A_matrix", or extract A.data (if A::IncidenceMatrix)
- evaluates BA with "calculate_BA_matrix", or extract BA.data (if A::BA_Matrix)
"""
function calculate_ABA_matrix(
    A::SparseArrays.SparseMatrixCSC{Int8, Int},
    BA::SparseArrays.SparseMatrixCSC{Float64, Int},
    ref_bus_positions::Set{Int})
    tmp = BA * A
    valid_ix = setdiff(1:size(tmp, 1), ref_bus_positions)
    return tmp[valid_ix, valid_ix]
end

"""
Return a sparse matrix given a dense one by dropping element whose absolute
value is above a certain tolerance.


# Arguments
- dense_array::Matrix{Float64}`:
        input matrix (e.g., PTDF matrix).
- `tol::Float64`:
        tolerance.
"""
function sparsify(dense_array::Matrix{Float64}, tol::Float64)
    m, n = size(dense_array)
    sparse_array = SparseArrays.spzeros(m, n)
    for j in 1:n, i in 1:m
        if abs(dense_array[i, j]) > tol
            sparse_array[i, j] = dense_array[i, j]
        end
    end
    return sparse_array
end

"""
Return a sparse vector given a dense one by dropping element whose absolute
value is above a certain tolerance.


# Arguments
- dense_array::Vector{Float64}`:
        input vector (e.g., PTDF row from VirtualPTDF).
- `tol::Float64`:
        tolerance.
"""
function sparsify(dense_array::Vector{Float64}, tol::Float64)
    m = length(dense_array)
    sparse_array = SparseArrays.spzeros(m)
    for i in 1:m
        if abs(dense_array[i]) > tol
            sparse_array[i] = dense_array[i]
        end
    end
    return sparse_array
end

"""
Takes the reference bus numbers and re-assigns the keys in the subnetwork dictionaries to use
the reference bus withing each subnetwork.
"""
function assign_reference_buses!(
    subnetworks::Dict{Int, Set{Int}},
    ref_buses::Vector{Int},
)
    if isempty(ref_buses) || length(ref_buses) != length(subnetworks)
        @warn "The reference bus positions are not consistent with the subnetworks. References buses will be assigned arbitrarily"
        return deepcopy(subnetworks)
    end
    bus_groups = Dict{Int, Set{Int}}()
    for (bus_key, subnetwork_buses) in subnetworks
        ref_bus = intersect(ref_buses, subnetwork_buses)
        if length(ref_bus) == 1
            bus_groups[first(ref_bus)] = pop!(subnetworks, bus_key)
        elseif length(ref_bus) == 0
            @warn "No reference bus in the subnetwork associated with bus $bus_key. Reference bus assigned arbitrarily"
        elseif length(ref_bus) > 1
            error(
                "More than one reference bus in the subnetwork associated with bus $bus_key",
            )
        else
            @assert false
        end
    end
    return bus_groups
end

function assign_reference_buses!(
    subnetworks::Dict{Int, Set{Int}},
    ref_bus_positions::Set{Int},
    bus_lookup::Dict{Int, Int},
)
    ref_buses = [k for (k, v) in bus_lookup if v in ref_bus_positions]
    return assign_reference_buses!(subnetworks, ref_buses)
end

"""
Finds the subnetworks present in the considered System. This is evaluated by taking
a the ABA or Adjacency Matrix.

# Arguments
- `M::SparseArrays.SparseMatrixCSC`:
        input sparse matrix.
- `bus_numbers::Vector{Int}`:
        vector containing the indices of the system's buses.
"""
function find_subnetworks(M::SparseArrays.SparseMatrixCSC, bus_numbers::Vector{Int})
    rows = SparseArrays.rowvals(M)
    touched = Set{Int}()
    subnetworks = Dict{Int, Set{Int}}()
    for (ix, bus_number) in enumerate(bus_numbers)
        neighbors = SparseArrays.nzrange(M, ix)
        if length(neighbors) <= 1
            @warn "Bus $bus_number is islanded"
            subnetworks[bus_number] = Set{Int}(bus_number)
            continue
        end
        for j in SparseArrays.nzrange(M, ix)
            row_ix = rows[j]
            if bus_number ∉ touched
                push!(touched, bus_number)
                subnetworks[bus_number] = Set{Int}(bus_number)
                _dfs(row_ix, M, bus_numbers, subnetworks[bus_number], touched)
            end
        end
    end
    return subnetworks
end

function _dfs(
    index::Int,
    M::SparseArrays.SparseMatrixCSC,
    bus_numbers::Vector{Int},
    bus_group::Set{Int},
    touched::Set{Int},
)
    rows = SparseArrays.rowvals(M)
    for j in SparseArrays.nzrange(M, index)
        row_ix = rows[j]
        if bus_numbers[row_ix] ∉ touched
            push!(touched, bus_numbers[row_ix])
            push!(bus_group, bus_numbers[row_ix])
            _dfs(row_ix, M, bus_numbers, bus_group, touched)
        end
    end
    return
end
