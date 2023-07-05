"""
Gets the AC branches from a given Systems.
"""
function get_ac_branches(sys::PSY.System)
    collection = Vector{PSY.ACBranch}()
    for br in PSY.get_components(PSY.get_available, PSY.ACBranch, sys)
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
        push!(collection, br)
    end
    return sort!(collection;
        by = x -> (PSY.get_number(PSY.get_arc(x).from), PSY.get_number(PSY.get_arc(x).to)),
    )
end

"""
Gets the non-isolated buses from a given System
"""
function get_buses(sys::PSY.System)::Vector{PSY.ACBus}
    return sort!(
        collect(
            PSY.get_components(x -> PSY.get_bustype(x) != ACBusTypes.ISOLATED, PSY.ACBus, sys),
        );
        by = x -> PSY.get_number(x),
    )
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
        (fr_b, to_b) = get_bus_indices(b, bus_lookup)
        a[fr_b, to_b] = 1
        a[to_b, fr_b] = -1
        a[fr_b, fr_b] = 1
        a[to_b, to_b] = 1
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
        to the refence buses
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
!!! MISSING DOCUMENTATION !!!
"""
function assing_reference_buses(
    subnetworks::Dict{Int, Set{Int}},
    ref_bus_positions::Set{Int},
)
    if isempty(ref_bus_positions) || length(ref_bus_positions) != length(subnetworks)
        @warn "The reference bus positions are not consistent with the subnetworks. References buses will be assigned arbitrarily"
        return deepcopy(subnetworks)
    end
    bus_groups = Dict{Int, Set{Int}}()
    for (bus_key, subnetwork_buses) in subnetworks
        ref_bus = intersect(ref_bus_positions, subnetwork_buses)
        if length(ref_bus) == 1
            bus_groups[first(ref_bus)] = pop!(subnetworks, bus_key)
            continue
        elseif length(ref_bus) == 0
            @warn "No reference bus in the subnetwork associated with bus $bus_key. References buses will be assigned arbitrarily"
            return subnetworks
        elseif length(ref_bus) > 1
            # TODO: still to implement
            error(
                "More than one reference bus in the subnetwork associated with bus $bus_key",
            )
        else
            @assert false
        end
    end
    return bus_groups
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
        if length(neighbors) < 1
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
