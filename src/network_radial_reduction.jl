struct RadialNetworkReduction
    bus_reduction_map::Dict{Int, Set{Int}}
    reverse_bus_search_map::Dict{Int, Int}
    radial_branches::Set{String}
    meshed_branches::Set{String}
end

get_bus_reduction_map(rb::RadialNetworkReduction) = rb.bus_reduction_map
get_reverse_bus_search_map(rb::RadialNetworkReduction) = rb.reverse_bus_search_map
get_radial_branches(rb::RadialNetworkReduction) = rb.radial_branches
get_meshed_branches(rb::RadialNetworkReduction) = rb.meshed_branches

function Base.isempty(rb::RadialNetworkReduction)
    if !isempty(rb.bus_reduction_map)
        return false
    end

    if !isempty(rb.bus_reduction_map)
        return false
    end

    if !isempty(rb.reverse_bus_search_map)
        return false
    end
    return true
end

function RadialNetworkReduction(;
    bus_reduction_map::Dict{Int, Set{Int}} = Dict{Int, Set{Int}}(),
    reverse_bus_search_map::Dict{Int, Int} = Dict{Int, Int}(),
    radial_branches::Set{String} = Set{String}(),
    meshed_branches::Set{String} = Set{String}())
    return RadialNetworkReduction(
        bus_reduction_map,
        reverse_bus_search_map,
        radial_branches,
        meshed_branches,
    )
end

"""
Structure to collect information the branches in the system that are radial and can be
ignored in the models.

# Arguments
- `A::IncidenceMatrix`: IncidenceMatrix

"""
function RadialNetworkReduction(A::IncidenceMatrix)
    return calculate_radial_branches(A.data, A.lookup[1], A.lookup[2], A.ref_bus_positions)
end

"""
Structure to collect information the branches in the system that are radial and can be
ignored in the models.

# Arguments
- `sys::PSY.System`: System Data
"""

function RadialNetworkReduction(sys::PSY.System)
    return RadialNetworkReduction(IncidenceMatrix(sys))
end

function _find_upstream_bus(
    A::SparseArrays.SparseMatrixCSC{Int8, Int64},
    j::Int,
    reverse_line_map::Dict{Int64, String},
    radial_branches::Set{String},
    reduced_buses::Set{Int},
    reverse_bus_map::Dict{Int, Int},
)
    row_ix = A.rowval[A.colptr[j]]
    parent = setdiff!(A[row_ix, :].nzind, j)[1]
    if reverse_bus_map[parent] ∈ reduced_buses
        row_ix = A.rowval[A.colptr[j] + 1]
        parent = setdiff!(A[row_ix, :].nzind, j)[1]
    end
    push!(radial_branches, reverse_line_map[row_ix])
    return parent
end

function _new_parent(
    A::SparseArrays.SparseMatrixCSC{Int8, Int64},
    parent::Int,
    bus_reduction_map_index::Dict{Int, Set{Int}},
    reverse_line_map::Dict{Int64, String},
    radial_branches::Set{String},
    reverse_bus_map::Dict{Int, Int},
)
    if length(SparseArrays.nzrange(A, parent)) == 2
        parent_bus_number = reverse_bus_map[parent]
        new_parent_val = _find_upstream_bus(
            A,
            parent,
            reverse_line_map,
            radial_branches,
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
            reverse_line_map,
            radial_branches,
            reverse_bus_map,
        )
    end
    return
end

function _reverse_search(
    A::SparseArrays.SparseMatrixCSC{Int8, Int64},
    j::Int,
    bus_reduction_map_index::Dict{Int, Set{Int}},
    reverse_line_map::Dict{Int64, String},
    radial_branches::Set{String},
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
        reverse_line_map,
        radial_branches,
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
        reverse_line_map,
        radial_branches,
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
function calculate_radial_branches(
    A::SparseArrays.SparseMatrixCSC{Int8, Int64},
    line_map::Dict{String, Int},
    bus_map::Dict{Int, Int},
    ref_bus_positions::Set{Int},
)
    lk = ReentrantLock()
    buscount = length(bus_map)
    radial_branches = Set{String}()
    reverse_line_map = Dict(reverse(kv) for kv in line_map)
    reverse_bus_map = Dict(reverse(kv) for kv in bus_map)
    bus_reduction_map_index = Dict{Int, Set{Int}}(k => Set{Int}() for k in keys(bus_map))
    Threads.@threads for j in 1:buscount
        if length(SparseArrays.nzrange(A, j)) == 1
            lock(lk) do
                _reverse_search(
                    A,
                    j,
                    bus_reduction_map_index,
                    reverse_line_map,
                    radial_branches,
                    reverse_bus_map,
                    ref_bus_positions,
                )
            end
        end
    end
    reverse_bus_search_map = _make_reverse_bus_search_map(bus_reduction_map_index, buscount)
    meshed_branches = Set{String}()
    for k in keys(line_map)
        if k in radial_branches
            continue
        end
        push!(meshed_branches, k)
    end
    return RadialNetworkReduction(
        bus_reduction_map_index,
        reverse_bus_search_map,
        radial_branches,
        meshed_branches,
    )
end

"""
Interface to obtain the parent bus number of a reduced bus when radial branches are eliminated

# Arguments
- `rb::RadialNetworkReduction`: RadialNetworkReduction object
- `bus_number::Int`: Bus number of the reduced bus
"""
function get_mapped_bus_number(rb::RadialNetworkReduction, bus_number::Int)
    if isempty(rb)
        return bus_number
    end
    return get(rb.reverse_bus_search_map, bus_number, bus_number)
end

"""
Interface to obtain the parent bus number of a reduced bus when radial branches are eliminated

# Arguments
- `rb::RadialNetworkReduction`: RadialNetworkReduction object
- `bus::ACBus`: Reduced bus
"""
function get_mapped_bus_number(rb::RadialNetworkReduction, bus::PSY.ACBus)
    return get_mapped_bus_number(rb, PSY.get_number(bus))
end

##############################################################################
########################### Auxiliary functions ##############################
##############################################################################

function isequal(
    rb1::RadialNetworkReduction,
    rb2::RadialNetworkReduction,
)
    for field in fieldnames(typeof(rb1))
        if getfield(rb1, field) != getfield(rb2, field)
            return false
        end
    end
    return true
end
