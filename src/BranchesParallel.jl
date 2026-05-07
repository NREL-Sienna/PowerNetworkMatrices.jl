mutable struct BranchesParallel{T <: PSY.ACTransmission} <: PSY.ACTransmission
    branches::Vector{T}
    equivalent_ybus::Union{Matrix{YBUS_ELTYPE}, Nothing}
end

function BranchesParallel(branches::Vector{T}) where {T <: PSY.ACTransmission}
    BranchesParallel(branches, nothing)
end
# Constructor for the mixed types
function BranchesParallel(branches::Vector{PSY.ACTransmission})
    return BranchesParallel{PSY.ACTransmission}(branches, nothing)
end

function add_branch!(bp::BranchesParallel{T}, branch::T) where {T <: PSY.ACTransmission}
    push!(bp.branches, branch)
end

get_branch_type(::BranchesParallel{T}) where {T <: PSY.ACTransmission} = T

function get_name(bp::BranchesParallel{T}) where {T <: PSY.ACTransmission}
    base_string = _longest_starting_substring(PSY.get_name.(bp.branches)...)
    if isempty(base_string)
        base_string = join(PSY.get_name.(bp.branches), "_") * "_"
    end
    return base_string *= "double_circuit"
end

function _longest_starting_substring(branch_names...)
    first_name = first(branch_names)
    n_chars = minimum(length.(branch_names))
    n_branches = length(branch_names)
    for ix in 1:n_chars
        for jx in 2:n_branches
            if branch_names[jx][ix] != first_name[ix]
                if ix == 1
                    return ""
                else
                    return first_name[1:(ix - 1)]
                end
            end
        end
    end
    return first_name[1:n_chars]
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
    return sum(get_series_susceptance(branch) for branch in segment.branches)
end

function get_equivalent_physical_branch_parameters(bp::BranchesParallel)
    if isnothing(bp.equivalent_ybus)
        populate_equivalent_ybus!(bp)
    end
    equivalent_ybus = bp.equivalent_ybus
    return _get_equivalent_physical_branch_parameters(equivalent_ybus)
end

function populate_equivalent_ybus!(bp::BranchesParallel)
    Y11, Y12, Y21, Y22 = ybus_branch_entries(bp)
    bp.equivalent_ybus = YBUS_ELTYPE[Y11 Y12; Y21 Y22]
    return
end

"""
    get_equivalent_rating(
        bp::BranchesParallel; 
        method::Symbol = :sum, 
        weighting::Symbol = :admittance_weighted
        )

Calculate the equivalent rating for a group of parallel branches.

Two orthogonal keyword arguments control the calculation:

## `weighting` — how flow is distributed across parallel branches

- `:admittance_weighted` (default): Each branch carries a fraction of the total flow
  proportional to its series susceptance ``f_i = b_i / \\sum_k b_k``, 
  where ``b_i`` is the series susceptance of branch i.
  This reflects the physical behavior of parallel circuits, where flow 
  distributes according to branch admittances.

- `:arithmetic`: All branches are treated as carrying equal fractions of total flow
  (uniform weighting).

## `method` — how individual branch limits are aggregated

- `:sum` (default): For `:admittance_weighted`, returns the total interface capacity
  limited by the first branch to reach its thermal limit (bottleneck formula):

  ``S_{\\max} = \\min_i \\left( \\frac{S_{\\text{limit},i}}{f_i} \\right)``

  For `:arithmetic`, returns the simple sum of individual ratings.

- `:average`: For `:admittance_weighted`, returns the susceptance-weighted average
  of individual ratings ``\\sum_i f_i \\cdot S_{\\text{limit},i}``.
  For `:arithmetic`, returns the arithmetic mean of individual ratings.

# Arguments
- `bp::BranchesParallel`: The parallel branch group.

# Keywords
- `method::Symbol = :sum`: Aggregation method. Valid values: `:sum`, `:average`.
- `weighting::Symbol = :admittance_weighted`: Flow weighting scheme. Valid values:
  `:admittance_weighted`, `:arithmetic`.
"""
function get_equivalent_rating(
    bp::BranchesParallel;
    method::Symbol = :sum,
    weighting::Symbol = :admittance_weighted,
)
    if weighting === :admittance_weighted
        if method === :sum
            # Total interface capacity limited by the first branch to reach its thermal limit.
            return minimum(
                get_equivalent_rating(br) /
                compute_parallel_multiplier(bp, PSY.get_name(br))
                for br in bp.branches
            )
        elseif method === :average
            # Susceptance-weighted average of individual ratings.
            return sum(
                compute_parallel_multiplier(bp, PSY.get_name(br)) *
                get_equivalent_rating(br)
                for br in bp.branches
            )
        else
            throw(
                ArgumentError(
                    "Unknown method: $(method). Valid options are :sum, :average.",
                ),
            )
        end
    elseif weighting === :arithmetic
        if method === :sum
            return sum(get_equivalent_rating(branch) for branch in bp.branches)
        elseif method === :average
            return sum(get_equivalent_rating(branch) for branch in bp.branches) /
                   length(bp.branches)
        else
            throw(
                ArgumentError(
                    "Unknown method: $(method). Valid options are :sum, :average.",
                ),
            )
        end
    else
        throw(
            ArgumentError(
                "Unknown weighting: $(weighting). Valid options are :admittance_weighted, :arithmetic.",
            ),
        )
    end
end

"""
    get_equivalent_emergency_rating(bp::BranchesParallel)

Calculate the total emergency rating for branches in parallel.
For parallel circuits, the emergency rating is the sum of individual emergency ratings divided by the number of circuits.
This provides a conservative estimate that accounts for potential overestimation of total capacity.
"""
function get_equivalent_emergency_rating(bp::BranchesParallel)
    equivalent_rating = 0.0
    for branch in bp.branches
        rating_b = get_equivalent_emergency_rating(branch)
        equivalent_rating += rating_b
    end
    return equivalent_rating # In Emergency conditions, we should consider the full capacity
end

"""
    get_equivalent_available(bp::BranchesParallel)

Get the availability status for parallel branches.
All branches in parallel must be available for the parallel circuit to be available.
"""
function get_equivalent_available(bp::BranchesParallel)
    # All branches must be available
    return all(PSY.get_available(branch) for branch in bp.branches)
end

"""
    get_equivalent_α(bp::BranchesParallel)

Get the phase angle shift for parallel branches.
Returns the average phase angle shift across all parallel branches.
Returns 0.0 if branches don't support phase angle shift (e.g., lines).
"""
function get_equivalent_α(bp::BranchesParallel)
    # Need to check the PS books
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

function add_to_map(
    double_circuit::BranchesParallel{T},
    filters::Dict,
) where {T <: PSY.ACTransmission}
    isempty(filters) && return true
    if isabstracttype(T)
        @warn "Parallel circuit contains mixed branch types, filters might be applied to more components than intended. Use Logging.Debug for additional information."
        @debug "Parallel circuit branch types: $(typeof.(double_circuit.branches))"
        @debug "Parallel circuit branch names: $(PSY.get_name.(double_circuit.branches))"
        for branch in double_circuit.branches
            filter = get(filters, typeof(branch), x -> true)
            if !filter(branch)
                return false
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

function Base.:(==)(a::BranchesParallel, b::BranchesParallel)
    return a.branches == b.branches
end

function Base.show(io::IO, x::MIME{Symbol("text/plain")}, y::BranchesParallel)
    show(io, x, y.branches)
end

is_a_reduction(::BranchesParallel) = true
