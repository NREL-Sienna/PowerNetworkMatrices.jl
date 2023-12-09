struct RadialBranches
    bus_reduction_map::Dict{Int, Set{Int}}
    radial_branches::Set{String}
end

get_bus_reduction_map(rb::RadialBranches) = rb.bus_reduction_map
get_radial_branches(rb::RadialBranches) = rb.radial_branches

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
        new_set = push!(pop!(bus_reduction_map_index, parent_bus_number), parent_bus_number)
        bus_reduction_map_index[new_parent_bus_number] = new_set
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
)
    j_bus_number = reverse_bus_map[j]
    reducion_set = Set{Int}(j_bus_number)
    parent = _find_upstream_bus(
        A,
        j,
        reverse_line_map,
        radial_branches,
        reducion_set,
        reverse_bus_map,
    )
    parent_bus_number = reverse_bus_map[parent]
    bus_reduction_map_index[parent_bus_number] = reducion_set
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

function RadialBranches(A::IncidenceMatrix)
    lk = ReentrantLock()
    buscount = size(A, 2)
    bus_reduction_map_index = Dict{Int, Set{Int}}()
    radial_branches = Set{String}()
    reverse_line_map = Dict(reverse(kv) for kv in A.lookup[1])
    reverse_bus_map = Dict(reverse(kv) for kv in A.lookup[2])
    Threads.@threads for j in 1:buscount
        if length(SparseArrays.nzrange(A.data, j)) == 1
            lock(lk) do
                _reverse_search(
                    A.data,
                    j,
                    bus_reduction_map_index,
                    reverse_line_map,
                    radial_branches,
                    reverse_bus_map,
                )
            end
        end
    end

    return RadialBranches(bus_reduction_map_index, radial_branches)
end
