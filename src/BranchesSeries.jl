mutable struct BranchesSeries <: PSY.ACTransmission
    branches::Dict{DataType, Vector{<:PSY.ACTransmission}}
    needs_insertion_order::Bool
    insertion_order::Vector{Tuple{DataType, Int}}
    segment_orientations::Vector{Symbol}
end

BranchesSeries() = BranchesSeries(
    Dict{DataType, Vector{<:PSY.ACTransmission}}(),
    false,
    Vector{Tuple{DataType, Int}}(),
    Vector{Symbol}(),
)

function add_branch!(
    bs::BranchesSeries,
    branch::T,
    orientation,
) where {T <: PSY.ACTransmission}
    push!(bs.segment_orientations, orientation)
    if isempty(bs.branches)
        # add the branch just once and return
        push!(get!(bs.branches, T, Vector{T}()), branch)
        return
    end

    if haskey(bs.branches, T) && !bs.needs_insertion_order
        push!(bs.branches[T], branch)
    elseif !haskey(bs.branches, T) && isempty(bs.insertion_order)
        bs.needs_insertion_order = true
        @assert length(keys(bs.branches)) == 1
        for (existing_type, existing_branches) in bs.branches
            for i in eachindex(existing_branches)
                push!(bs.insertion_order, (existing_type, i))
            end
        end
        push!(get!(bs.branches, T, Vector{T}()), branch)
        push!(bs.insertion_order, (T, 1))
    elseif !haskey(bs.branches, T) && !isempty(bs.insertion_order)
        push!(get!(bs.branches, T, Vector{T}()), branch)
        push!(bs.insertion_order, (T, length(bs.branches[T])))
    else
        push!(bs.branches[T], branch)
        push!(bs.insertion_order, (T, length(bs.branches[T])))
    end
    return
end

# Iteration support
function Base.iterate(bs::BranchesSeries)
    if isempty(bs.branches)
        return nothing
    end

    # Single key case - iterate over the single vector directly
    if bs.needs_insertion_order
        # Multi-key case - use insertion_order
        if isempty(bs.insertion_order)
            return nothing
        end

        type, idx = bs.insertion_order[1]
        branch = bs.branches[type][idx]
        return (branch, (1, nothing))
    else
        single_vector = first(values(bs.branches))
        if isempty(single_vector)
            return nothing
        end
        return (single_vector[1], (1, single_vector))
    end
end

function Base.iterate(bs::BranchesSeries, state)
    position, vector_cache = state

    if bs.needs_insertion_order
        # Multi-key iteration using insertion_order
        next_position = position + 1
        if next_position > length(bs.insertion_order)
            return nothing
        end
        type, idx = bs.insertion_order[next_position]
        branch = bs.branches[type][idx]
        return (branch, (next_position, nothing))
    else
        # Single key iteration
        next_idx = position + 1
        if next_idx > length(vector_cache)
            return nothing
        end
        return (vector_cache[next_idx], (next_idx, vector_cache))
    end
end

Base.length(bs::BranchesSeries) =
    if bs.needs_insertion_order
        length(bs.insertion_order)
    else
        sum(length(v) for v in values(bs.branches))
    end

Base.eltype(::Type{BranchesSeries}) = PSY.ACTransmission

function get_series_susceptance(series_chain::BranchesSeries)
    series_susceptances_sum = sum(inv(get_series_susceptance(x)) for x in series_chain)
    total_susceptance = 1 / series_susceptances_sum
    return total_susceptance
end

"""
    get_equivalent_physical_branch_parameters(bs::BranchesSeries)

Calculate the physical parameters for branches in series.
This method computes an equivalent Ybus matrix for the, reduces the internal nodes, 
    and then computes equivalent physical parameters for the equivalent ybus. 
"""
function get_equivalent_physical_branch_parameters(bs::BranchesSeries)
    ybus_isolated_d2_chain = _build_chain_ybus(bs)
    ybus_boundary_isolated_d2_chain = _reduce_internal_nodes(ybus_isolated_d2_chain)
    y_11, y_12, y_21, y_22 = ybus_boundary_isolated_d2_chain
    if isapprox(y_12, y_21)
        tap = 1.0
        shift = 0.0
    else
        error("Non-symmetric components not yet implemented for BranchesSeries")
    end
    y_l = y_12 * -1
    z_12 = 1 / y_l
    r = real(z_12)
    x = imag(z_12)
    g_from = real(y_11 - y_l)
    b_from = imag(y_11 - y_l)
    g_to = real(y_22 - y_l)
    b_to = imag(y_22 - y_l)
    return r, x, g_from, b_from, g_to, b_to, tap, shift
end

"""
    get_equivalent_rating(bs::BranchesSeries)

Calculate the rating for branches in series.
For series circuits, the rating is limited by the weakest link: Rating_total = min(Rating1, Rating2, ..., Ratingn)
"""
function get_equivalent_rating(bs::BranchesSeries)
    # Minimum rating for series branches (weakest link)
    return minimum(PSY.get_rating(branch) for branch in bs)
end

"""
    get_equivalent_available(bs::BranchesSeries)

Get the availability status for series branches.
All branches in series must be available for the series circuit to be available.
"""
function get_equivalent_available(bs::BranchesSeries)
    # All branches must be available
    return all(PSY.get_available(branch) for branch in bs)
end

"""
    get_equivalent_α(bs::BranchesSeries)

Get the phase angle shift for series branches.
Returns the sum of phase angle shifts across all series branches.
Returns 0.0 if branches don't support phase angle shift (e.g., lines).
"""
function get_equivalent_α(bs::BranchesSeries)
    # Need to check how to develop this one
end

function add_to_map(series_circuit::BranchesSeries, filters::Dict)
    if isempty(filters)
        return true
    end

    if series_circuit.needs_insertion_order
        if isempty(intersect(keys(series_circuit.branches), keys(filters)))
            return true
        end

        @warn "Series circuit contains mixed branch types, filters might be applied to more components than intended. Use Logging.Debug for additional information."
        @debug "Series circuit branch types: $(keys(series_circuit.branches))"
        @debug "Series circuit branch names: $(vcat([PSY.get_name.(v) for (k , v) in series_circuit.branches]))"
        for (branch_type, branch_list) in series_circuit.branches
            filter = get(filters, branch_type, x -> true)
            for device in branch_list
                if !filter(device)
                    return false
                end
            end
        end
        return true
    else
        filter = filters[first(keys(series_circuit.branches))]
        return all([filter(device) for device in first(values(series_circuit.branches))])
    end
    error("Invalid condition reached in add_to_map for BranchesSeries")
end

function Base.:(==)(a::BranchesSeries, b::BranchesSeries)
    return a.branches == b.branches
end

function Base.show(io::IO, x::MIME{Symbol("text/plain")}, y::BranchesSeries)
    show(io, x, y.branches)
end
