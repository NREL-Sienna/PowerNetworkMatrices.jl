struct BranchesParallel{T <: PSY.ACTransmission} <: PSY.ACTransmission
    branches::Set{T}
end

function BranchesParallel(branches::Vector{T}) where {T <: PSY.ACTransmission}
    return BranchesParallel{T}(Set{T}(branches))
end

function BranchesParallel(branches::Vector)
    return BranchesParallel{ACTransmission}(Set{ACTransmission}(branches))
end

function add_branch!(bp::BranchesParallel{T}, branch::T) where {T <: PSY.ACTransmission}
    push!(bp.branches, branch)
end

get_branch_type(::BranchesParallel{T}) where {T <: PSY.ACTransmission} = T

function get_name(bp::BranchesParallel{T}) where {T <: PSY.ACTransmission}
    base_string = join(intersect(PSY.get_name.(bp.branches)...))
    return base_string *= "double_circuit"
end

function compute_parallel_multiplier(
    parallel_branch_set::BranchesParallel,
    branch_name::String,
)
    b_total = 0.0
    b_branch = 0.0
    for br in parallel_branch_set
        if PSY.get_name(br) == branch_name
            b_branch += PSY.get_series_susceptance(br)
        end
        b_total += PSY.get_series_susceptance(br)
    end
    return b_branch / b_total
end

function get_series_susceptance(segment::BranchesParallel)
    return sum([get_series_susceptance(branch) for branch in segment.branches])
end

function Base.iterate(bp::BranchesParallel)
    return iterate(bp.branches)
end

function Base.iterate(bp::BranchesParallel, state)
    return iterate(bp.branches, state)
end

function Base.length(bp::BranchesParallel)
    return length(bp.branches)
end

function add_to_map(double_circuit::BranchesParallel{T}, filters::Dict) where {T <: PSY.ACTransmission}
    if isabstracttype(T)
        @warn "Parallel circuit contains mixed branch types, filters might be applied to more components than intended. Use Logging.Debug for additional information."
        @debug "Parallel circuit branch types: $(keys(double_circuit.branches))"
        @debug "Parallel circuit branch names: $(vcat([PSY.get_name.(v) for (k , v) in double_circuit.branches]))"
        for branch in double_circuit.branches
            filter = get(filters, typeof(branch), x -> true)
            for device in branch
                if !filter(device)
                    return false
                end
            end
        end
        return true
    end

    if !haskey(filters, T)
        return true
    else
        return any([filters[T](device) for device in double_circuit])
    end
    error("Invalid condition reached in add_to_map for BranchesParallel")
end
