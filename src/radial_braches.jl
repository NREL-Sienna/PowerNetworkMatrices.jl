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
