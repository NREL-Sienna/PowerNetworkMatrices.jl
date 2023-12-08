struct RadialBranches
    bus_reduction_map_index::Dict{Int, Set{Int}}
    radial_branches::Vector{String}
end

function RadialBranches(A::IncidenceMatrix)
    buscount = size(A, 2)
    bus_reduction_map_index = Dict{Int, Set{Int}}()
    radial_branches = Set{String}()
    reverse_line_map = Dict(reverse(kv) for kv in A.lookup[1])
    for j in 1:buscount
        @show length(SparseArrays.nzrange(A.data, j))
        if length(SparseArrays.nzrange(A.data, j)) == 1
            push!(radial_branches, reverse_line_map[A.data.rowval[j]])

        end
    end

    return RadialBranches(bus_reduction_map_index, radial_branches)
end

function _find_upstream_bus(A::SparseMatrixCSC{Int8, Int64}, j::Int, reverse_line_map::Dict{Int64, String}, radial_branches::Set{String}, reduced_buses::Set{Int})
    row_ix = A.rowval[A.colptr[j]]
    parent = setdiff!(A[row_ix, :].nzind, j)[1]
    if parent ∈ reduced_buses
        row_ix = A.rowval[A.colptr[j] + 1]
        parent = setdiff!(A[row_ix, :].nzind, j)[1]
    end
    push!(radial_branches, reverse_line_map[row_ix])
    return parent
end

function _new_parent(A::SparseMatrixCSC{Int8, Int64}, parent::Int, bus_reduction_map_index::Dict{Int, Set{Int}}, reverse_line_map::Dict{Int64, String}, radial_branches::Set{String})
    if length(SparseArrays.nzrange(A, parent)) == 2
        new_parent_val = find_upstream_bus(A, parent, reverse_line_map, radial_branches, bus_reduction_map_index[parent])
        new_set = push!(pop!(bus_reduction_map_index, parent), parent)
        bus_reduction_map_index[new_parent_val] = new_set
        _new_parent(A, new_parent_val, bus_reduction_map_index, reverse_line_map, radial_branches)
    end
    return
end

function _reverse_search(A::SparseMatrixCSC{Int8, Int64}, j::Int, bus_reduction_map_index::Dict{Int, Set{Int}}, reverse_line_map, radial_branches::Set{String})
    reducion_set = Set{Int}(j)
    parent = find_upstream_bus(A, j, reverse_line_map, radial_branches, reducion_set)
    bus_reduction_map_index[parent] = reducion_set
    _new_parent(A, parent, bus_reduction_map_index, reverse_line_map, radial_branches)
    return
end

function RadialBranches(A::IncidenceMatrix)
    lk = ReentrantLock()
    buscount = size(A, 2)
    bus_reduction_map_index = Dict{Int, Set{Int}}()
    radial_branches = Set{String}()
    reverse_line_map = Dict(reverse(kv) for kv in A.lookup[1])
    Threads.@threads for j in 1:buscount
        if length(SparseArrays.nzrange(A.data, j)) == 1
            lock(lk) do
            _reverse_search(A.data, j, bus_reduction_map_index, reverse_line_map, radial_branches)
        end
    end
    end

    return bus_reduction_map_index, radial_branches
end
