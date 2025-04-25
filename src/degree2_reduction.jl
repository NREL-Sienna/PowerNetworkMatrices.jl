function degree2_reduction(sys::PSY.System)

    irreducible_buses = Int[]
    for c in PSY.get_components(PSY.StaticInjection, sys)
        push!(irreducible_buses, PSY.get_number(PSY.get_bus(c)))
    end

    irreducible_indices = Set([p.lookup[1][i] for i in irreducible_buses])

    A = IncidenceMatrix(sys)

    I = Adjacency(sys)
    return degree2_reduction(irreducible_indices, A, I)
end

function degree2_reduction(irreducible_buses::Vector{Int}, A::Adjacency)

    degree_two_indexes = Int[]

    for i in 1:size(A.data, 1)
        nzrange = SparseArrays.nzrange(A.data, i)
        if length(SparseArrays.nzrange(p.data, i)) == 3
            push!(degree_two_indexes, i)
            @show rowvals(p.data)[nzrange]
            @show left_n = first(rowvals(p.data)[nzrange])
            @show right_n = last(rowvals(p.data)[nzrange])
            @show length(SparseArrays.nzrange(p.data, left_n))
            @show length(SparseArrays.nzrange(p.data, right_n))
        end
    end


    reducible_2d_buses_indexes = setdiff(degree_two_indexes, irreducible_indices)

end

function get_branch_joins(I::IncidenceMatrix,reducible_2d_buses_indexes::Vector{Int})
    branch_joins = []
    for i in reducible_2d_buses_indexes
        @show rowvals(A.data)[SparseArrays.nzrange(I.data, i)]
        @assert length(rowvals(A.data)[SparseArrays.nzrange(I.data, i)]) == 2 i
        push!(branch_joins, rowvals(A.data)[SparseArrays.nzrange(I.data, i)])
    end

    return branch_joins
end
