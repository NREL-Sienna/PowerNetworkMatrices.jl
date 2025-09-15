function _add_to_collection!(
    collection_br::Vector{PSY.ACTransmission},
    branch::PSY.ACTransmission,
)
    push!(collection_br, branch)
    return
end

function _add_to_collection!(
    collection_tr3w::Vector{PSY.ThreeWindingTransformer},
    transformer_tr3w::PSY.ThreeWindingTransformer,
)
    push!(collection_tr3w, transformer_tr3w)
    return
end

function get_bus_index(bus_no::Int, bus_lookup::Dict{Int, Int}, nr::NetworkReductionData)
    if haskey(nr.reverse_bus_search_map, bus_no)
        return bus_lookup[nr.reverse_bus_search_map[bus_no]]
    else
        return bus_lookup[bus_no]
    end
end

function get_bus_index(
    dev::PSY.Component,
    bus_lookup::Dict{Int, Int},
    nr::NetworkReductionData,
)
    bus_number = PSY.get_number(PSY.get_bus(dev))
    return get_bus_index(bus_number, bus_lookup, nr)
end

function get_bus_indices(
    arc::PSY.Arc,
    bus_lookup::Dict{Int, Int},
    nr::NetworkReductionData,
)
    check_arc_validity(arc, IS.get_name(arc))
    reverse_bus_search_map = get_reverse_bus_search_map(nr)
    fr_bus_number = PSY.get_number(PSY.get_from(arc))
    if haskey(reverse_bus_search_map, fr_bus_number)
        fr_bus_number_reduced = reverse_bus_search_map[fr_bus_number]
    else
        fr_bus_number_reduced = fr_bus_number
    end
    fr_bus_ix = bus_lookup[fr_bus_number_reduced]

    to_bus_number = PSY.get_number(PSY.get_to(arc))
    if haskey(reverse_bus_search_map, to_bus_number)
        to_bus_number_reduced = reverse_bus_search_map[to_bus_number]
    else
        to_bus_number_reduced = to_bus_number
    end
    to_bus_ix = bus_lookup[to_bus_number_reduced]
    return fr_bus_ix, to_bus_ix
end

function check_arc_validity(arc::PSY.Arc, name::String)
    if PSY.get_bustype(PSY.get_from(arc)) == ACBusTypes.ISOLATED
        throw(
            IS.ConflictingInputsError(
                "Branch or arc $(name) is set available and connected to isolated bus " *
                "$(IS.get_name(PSY.get_from(arc)))",
            ),
        )
    end
    if PSY.get_bustype(PSY.get_to(arc)) == ACBusTypes.ISOLATED
        throw(
            IS.ConflictingInputsError(
                "Branch or arc $(name) is set available and connected to isolated bus " *
                "$(IS.get_name(PSY.get_to(arc)))",
            ),
        )
    end
    return
end

function get_arc_tuple(arc::PSY.Arc, nr::NetworkReductionData)
    reverse_bus_search_map = get_reverse_bus_search_map(nr)
    arc_tuple_original = get_arc_tuple(arc)
    return (
        get(reverse_bus_search_map, arc_tuple_original[1], arc_tuple_original[1]),
        get(reverse_bus_search_map, arc_tuple_original[2], arc_tuple_original[2]),
    )
end
function get_arc_tuple(br::PSY.ACTransmission, nr::NetworkReductionData)
    get_arc_tuple(PSY.get_arc(br), nr)
end

#Assumes set of branches are in parallel 
function get_arc_tuple(br::Set{PSY.ACTransmission}, nr::NetworkReductionData)
    get_arc_tuple(PSY.get_arc(first(br)), nr)
end

function get_arc_tuple(br::Set{PSY.ACTransmission})
    return get_arc_tuple(PSY.get_arc(first(br)))
end

function get_arc_tuple(br::PSY.ACTransmission)
    return get_arc_tuple(PSY.get_arc(br))
end

get_arc_tuple(arc::PSY.Arc) =
    (PSY.get_number(PSY.get_from(arc)), PSY.get_number(PSY.get_to(arc)))

function get_arc_tuple(tr3W_tuple::Tuple{T, Int}) where {T <: PSY.ThreeWindingTransformer}
    t3W, arc_number = tr3W_tuple
    if arc_number == 1
        return (
            PSY.get_number(PSY.get_from(PSY.get_primary_star_arc(t3W))),
            PSY.get_number(PSY.get_to(PSY.get_primary_star_arc(t3W))),
        )
    elseif arc_number == 2
        return (
            PSY.get_number(PSY.get_from(PSY.get_secondary_star_arc(t3W))),
            PSY.get_number(PSY.get_to(PSY.get_secondary_star_arc(t3W))),
        )
    elseif arc_number == 3
        return (
            PSY.get_number(PSY.get_from(PSY.get_tertiary_star_arc(t3W))),
            PSY.get_number(PSY.get_to(PSY.get_tertiary_star_arc(t3W))),
        )
    end
end

function get_arc_tuple(
    tr3W_tuple::Tuple{T, Int},
    nr::NetworkReductionData,
) where {T <: PSY.ThreeWindingTransformer}
    reverse_bus_search_map = get_reverse_bus_search_map(nr)
    arc_tuple_original = get_arc_tuple(tr3W_tuple)
    return (
        get(reverse_bus_search_map, arc_tuple_original[1], arc_tuple_original[1]),
        get(reverse_bus_search_map, arc_tuple_original[2], arc_tuple_original[2]),
    )
end

"""
Gets the AC branches for a system
"""
function get_ac_branches(
    sys::PSY.System,
    radial_branches::Set{String} = Set{String}(),
)::Vector{PSY.ACTransmission}
    collection_br = Vector{PSY.ACTransmission}()
    for br in PSY.get_components(
        x -> PSY.get_available(x) && !(typeof(x) <: PSY.Transformer3W),
        PSY.ACTransmission,
        sys,
    )
        arc = PSY.get_arc(br)
        name = PSY.get_name(br)
        check_arc_validity(arc, name)

        if PSY.get_name(br) ∉ radial_branches
            _add_to_collection!(collection_br, br)
        end
    end

    return sort!(collection_br;
        by = get_arc_tuple,
    )
end

function get_switched_admittances(sys::PSY.System, reverse_bus_search_map)
    collection = Vector{PSY.SwitchedAdmittance}()
    for sa in
        collect(PSY.get_components(x -> PSY.get_available(x), PSY.SwitchedAdmittance, sys))
        if !haskey(reverse_bus_search_map, PSY.get_number(PSY.get_bus(sa)))
            push!(collection, sa)
        end
    end
    return collection
end

function get_fixed_admittances(sys::PSY.System, reverse_bus_search_map)
    collection = Vector{PSY.FixedAdmittance}()
    for sa in
        collect(PSY.get_components(x -> PSY.get_available(x), PSY.FixedAdmittance, sys))
        if !haskey(reverse_bus_search_map, PSY.get_number(PSY.get_bus(sa)))
            push!(collection, sa)
        end
    end
    return collection
end

function _add_branch_to_lookup!(
    branch_lookup::Dict{String, Int},
    ::Dict{String, Vector{String}},
    branch_type::Vector{DataType},
    branch::PSY.ACTransmission,
    branch_number::Int,
)
    branch_lookup[PSY.get_name(branch)] = branch_number
    push!(branch_type, typeof(branch))
    return
end

function _add_branch_to_lookup!(
    branch_lookup::Dict{String, Int},
    transformer_3w_lookup::Dict{String, Vector{String}},
    branch_type::Vector{DataType},
    branch::PSY.ThreeWindingTransformer,
    branch_number::Int,
)
    tr3w_name = PSY.get_name(branch)
    transformer_3w_lookup[tr3w_name] = Vector{String}(undef, 3)
    for (i, side) in enumerate(["primary", "secondary", "tertiary"])
        side_name = "$(tr3w_name)__$side"
        branch_lookup[side_name] = branch_number - 3 + i
        transformer_3w_lookup[tr3w_name][i] = side_name
        push!(branch_type, typeof(branch))
    end
    return
end

"""
Gets the non-isolated buses from a given System
"""
function get_buses(
    sys::PSY.System,
    bus_reduction_map::Dict{Int, Set{Int}} = Dict{Int, Set{Int}}(),
)::Vector{PSY.ACBus}
    leaf_buses = Set{PSY.Int}()
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
Validates that the user bus input is consistent with the ybus axes and the prior reductions.
Is used to check `irreducible_buses` for `Radial` and `DegreeTwo` reductions and `study_buses` for `WardReduction`.
"""
function validate_buses(A::AdjacencyMatrix, buses::Vector{Int})
    reverse_bus_search_map = A.network_reduction_data.reverse_bus_search_map
    for bus_no in buses
        reduced_bus_no = get(reverse_bus_search_map, bus_no, bus_no)
        if reduced_bus_no ∉ get_bus_axis(A)
            if bus_no == reduced_bus_no
                error(
                    "Invalid bus entry found: Bus $bus_no. Check your input data; this bus was not found in the admittance matrix.",
                )
            else
                error(
                    "Invalid bus entry found: Bus $bus_no. Check your input data; this bus was mapped to bus $reduced_bus_no in a prior reductions and not found in the admittance matrix.",
                )
            end
        end
    end
    return
end

"""
Convert the user input for irreducible_buses to a set of indices based on the Ybus lookup and the prior reductions.
"""
function get_irreducible_indices(A::AdjacencyMatrix, irreducible_buses::Vector{Int})
    reverse_bus_search_map = A.network_reduction_data.reverse_bus_search_map
    irreducible_indices = zeros(Int, length(irreducible_buses))
    for (ix, bus_no) in enumerate(irreducible_buses)
        reduced_bus_no = get(reverse_bus_search_map, bus_no, bus_no)
        irreducible_indices[ix] = A.lookup[1][reduced_bus_no]
    end
    return irreducible_indices
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
