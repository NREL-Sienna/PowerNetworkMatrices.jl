function get_reduction(
    A::IncidenceMatrix,
    sys::PSY.System,
    ::Val{NetworkReductionTypes.RADIAL},
)
    return get_radial_reduction(A)
end

"""
Builds a NetworkReduction by removing radially connected buses. 

# Arguments
- `sys::System`
"""
function get_radial_reduction(
    sys::PSY.System;
    prior_reduction::NetworkReduction = NetworkReduction(),
    exempt_buses::Vector{Int64} = Int64[],
)
    validate_reduction_type(
        NetworkReductionTypes.RADIAL,
        get_reduction_type(prior_reduction),
    )
    return get_radial_reduction(
        IncidenceMatrix(sys);
        prior_reduction = prior_reduction,
        exempt_buses = exempt_buses,
    )
end

"""
Builds a NetworkReduction by removing radially connected buses. 

# Arguments
- `A::IncidenceMatrix`: IncidenceMatrix
"""
function get_radial_reduction(
    A::IncidenceMatrix;
    exempt_buses::Vector{Int64} = Int64[],
)
    exempt_bus_numbers = union(A.ref_bus_numbers, Set(exempt_buses))
    exempt_bus_positions = Set([A.lookup[2][x] for x in exempt_bus_numbers])
    return calculate_radial_arcs(
        A.data,
        A.lookup[1],
        A.lookup[2],
        exempt_bus_positions,
    )
end

function _find_upstream_bus(
    A::SparseArrays.SparseMatrixCSC{Int8, Int64},
    j::Int,
    reverse_arc_map::Dict{Int64, Tuple{Int, Int}},
    radial_arcs::Set{Tuple{Int, Int}},
    reduced_buses::Set{Int},
    reverse_bus_map::Dict{Int, Int},
)
    row_ix = A.rowval[A.colptr[j]]
    parent = setdiff!(A[row_ix, :].nzind, j)[1]
    if reverse_bus_map[parent] ∈ reduced_buses
        row_ix = A.rowval[A.colptr[j] + 1]
        parent = setdiff!(A[row_ix, :].nzind, j)[1]
    end
    push!(radial_arcs, reverse_arc_map[row_ix])
    return parent
end

function _new_parent(
    A::SparseArrays.SparseMatrixCSC{Int8, Int64},
    parent::Int,
    bus_reduction_map_index::Dict{Int, Set{Int}},
    reverse_arc_map::Dict{Int64, Tuple{Int, Int}},
    radial_arcs::Set{Tuple{Int, Int}},
    reverse_bus_map::Dict{Int, Int},
)
    if length(SparseArrays.nzrange(A, parent)) == 2
        parent_bus_number = reverse_bus_map[parent]
        new_parent_val = _find_upstream_bus(
            A,
            parent,
            reverse_arc_map,
            radial_arcs,
            bus_reduction_map_index[parent_bus_number],
            reverse_bus_map,
        )
        new_parent_bus_number = reverse_bus_map[new_parent_val]
        # This check is meant to capture cases of a full radial network which can happen in
        # system with small islands that represent larger interconnected areas.
        if length(SparseArrays.nzrange(A, new_parent_val)) < 2
            @warn "Bus $parent_bus_number Parent $new_parent_bus_number is a leaf node, indicating there is an island."
            push!(bus_reduction_map_index[parent_bus_number], new_parent_bus_number)
            return
        end
        new_set = push!(pop!(bus_reduction_map_index, parent_bus_number), parent_bus_number)
        union!(bus_reduction_map_index[new_parent_bus_number], new_set)
        _new_parent(
            A,
            new_parent_val,
            bus_reduction_map_index,
            reverse_arc_map,
            radial_arcs,
            reverse_bus_map,
        )
    end
    return
end

function _reverse_search(
    A::SparseArrays.SparseMatrixCSC{Int8, Int64},
    j::Int,
    bus_reduction_map_index::Dict{Int, Set{Int}},
    reverse_arc_map::Dict{Int64, Tuple{Int, Int}},
    radial_arcs::Set{Tuple{Int, Int}},
    reverse_bus_map::Dict{Int, Int},
    ref_bus_positions::Set{Int},
)
    if j ∈ ref_bus_positions
        return
    end
    j_bus_number = reverse_bus_map[j]
    pop!(bus_reduction_map_index, j_bus_number)
    reduction_set = Set{Int}(j_bus_number)
    parent = _find_upstream_bus(
        A,
        j,
        reverse_arc_map,
        radial_arcs,
        reduction_set,
        reverse_bus_map,
    )
    parent_bus_number = reverse_bus_map[parent]
    union!(bus_reduction_map_index[parent_bus_number], reduction_set)
    if parent ∈ ref_bus_positions
        return
    end
    _new_parent(
        A,
        parent,
        bus_reduction_map_index,
        reverse_arc_map,
        radial_arcs,
        reverse_bus_map,
    )
    return
end

function _make_reverse_bus_search_map(bus_reduction_map::Dict{Int, Set{Int}}, n_buses::Int)
    map = Dict{Int, Int}()
    sizehint!(map, n_buses)
    for (parent, children) in bus_reduction_map
        for bus in children
            map[bus] = parent
        end
    end
    return map
end

"""
used to calculate the branches in the system that are radial and can be
ignored in the models by exploring the structure of the incidence matrix

# Arguments
- `A::SparseArrays.SparseMatrixCSC{Int8, Int64}`: Data from the IncidenceMatrix
- `line_map::Dict{String, Int}`: Map of Line Name to Matrix Index
- `bus_map::Dict{Int, Int}`: Map of Bus Name to Matrix Index
- `ref_bus_positions::Set{Int}`:
        Set containing the indexes of the columns of the BA matrix corresponding
        to the reference buses
"""
function calculate_radial_arcs(
    A::SparseArrays.SparseMatrixCSC{Int8, Int64},
    arc_map::Dict{Tuple{Int, Int}, Int},
    bus_map::Dict{Int, Int},
    ref_bus_positions::Set{Int},
)
    lk = ReentrantLock()
    buscount = length(bus_map)
    radial_arcs = Set{Tuple{Int, Int}}()
    reverse_arc_map = Dict(reverse(kv) for kv in arc_map)
    reverse_bus_map = Dict(reverse(kv) for kv in bus_map)
    bus_reduction_map_index = Dict{Int, Set{Int}}(k => Set{Int}() for k in keys(bus_map))
    Threads.@threads for j in 1:buscount
        if length(SparseArrays.nzrange(A, j)) == 1
            lock(lk) do
                _reverse_search(
                    A,
                    j,
                    bus_reduction_map_index,
                    reverse_arc_map,
                    radial_arcs,
                    reverse_bus_map,
                    ref_bus_positions,
                )
            end
        end
    end
    reverse_bus_search_map = _make_reverse_bus_search_map(bus_reduction_map_index, buscount)

    return NetworkReduction(
        bus_reduction_map_index,
        reverse_bus_search_map,
        Dict{Tuple{Int, Int}, PSY.Branch}(),
        Dict{PSY.Branch, Tuple{Int, Int}}(),
        Dict{Tuple{Int, Int}, Set{PSY.Branch}}(),
        Dict{PSY.Branch, Tuple{Int, Int}}(),
        Dict{Tuple{Int, Int}, Set{PSY.Branch}}(),
        Dict{PSY.Branch, Tuple{Int, Int}}(),
        Dict{Tuple{Int, Int}, Tuple{PSY.Transformer3W, Int}}(),
        Dict{Tuple{PSY.Transformer3W, Int}, Tuple{Int, Int}}(),
        Set{Int}(),
        radial_arcs,
        [NetworkReductionTypes.RADIAL],
    )
end

"""
Interface to obtain the parent bus number of a reduced bus when radial branches are eliminated

# Arguments
- `rb::NetworkReduction`: NetworkReduction object
- `bus_number::Int`: Bus number of the reduced bus
"""
function get_mapped_bus_number(rb::NetworkReduction, bus_number::Int)
    if isempty(rb)
        return bus_number
    end
    return get(rb.reverse_bus_search_map, bus_number, bus_number)
end

"""
Interface to obtain the parent bus number of a reduced bus when radial branches are eliminated

# Arguments
- `rb::NetworkReduction`: NetworkReduction object
- `bus::ACBus`: Reduced bus
"""
function get_mapped_bus_number(rb::NetworkReduction, bus::PSY.ACBus)
    return get_mapped_bus_number(rb, PSY.get_number(bus))
end
