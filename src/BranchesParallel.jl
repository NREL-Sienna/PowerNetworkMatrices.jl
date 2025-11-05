struct BranchesParallel{T <: PSY.ACTransmission} <: PSY.ACTransmission
    branches::Vector{T}
end

# Constructor for the mixed types
function BranchesParallel(branches::Vector{PSY.ACTransmission})
    return BranchesParallel{PSY.ACTransmission}(Vector{PSY.ACTransmission}(branches))
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

"""
    get_equivalent_r(bp::BranchesParallel)

Calculate the equivalent resistance for branches in parallel.
For parallel circuits, the equivalent impedance is: 1/Z_eq = 1/Z1 + 1/Z2 + ... + 1/Zn
where Z = R + jX. The equivalent resistance is the real part of Z_eq.
"""
function get_equivalent_r(bp::BranchesParallel)
    # Calculate equivalent impedance: 1/Z_eq = sum(1/Z_i)
    inv_z_sum =
        sum(inv(PSY.get_r(branch) + 1im * PSY.get_x(branch)) for branch in bp.branches)
    z_eq = inv(inv_z_sum)
    return real(z_eq)
end

"""
    get_equivalent_x(bp::BranchesParallel)

Calculate the equivalent reactance for branches in parallel.
For parallel circuits, the equivalent impedance is: 1/Z_eq = 1/Z1 + 1/Z2 + ... + 1/Zn
where Z = R + jX. The equivalent reactance is the imaginary part of Z_eq.
"""
function get_equivalent_x(bp::BranchesParallel)
    # Calculate equivalent impedance: 1/Z_eq = sum(1/Z_i)
    inv_z_sum =
        sum(inv(PSY.get_r(branch) + 1im * PSY.get_x(branch)) for branch in bp.branches)
    z_eq = inv(inv_z_sum)
    return imag(z_eq)
end

"""
    get_equivalent_b(bp::BranchesParallel)

Calculate the equivalent susceptance for branches in parallel.
For parallel circuits: B_total = (from = B1_from + B2_from + ..., to = B1_to + B2_to + ...)
Returns a NamedTuple with :from and :to fields.
"""
function get_equivalent_b(bp::BranchesParallel)
    # Direct sum for parallel susceptances
    b_from = sum(PSY.get_b(branch).from for branch in bp.branches)
    b_to = sum(PSY.get_b(branch).to for branch in bp.branches)
    return (from = b_from, to = b_to)
end

"""
    get_equivalent_g(bp::BranchesParallel)

Calculate the equivalent conductance for branches in parallel.
For parallel circuits: G_total = (from = G1_from + G2_from + ..., to = G1_to + G2_to + ...)
Returns a NamedTuple with :from and :to fields.
"""
function get_equivalent_g(bp::BranchesParallel)
    # Direct sum for parallel conductances
    g_from = sum(PSY.get_g(branch).from for branch in bp.branches)
    g_to = sum(PSY.get_g(branch).to for branch in bp.branches)
    return (from = g_from, to = g_to)
end

"""
    get_equivalent_rating(bp::BranchesParallel)

Calculate the total rating for branches in parallel.
For parallel circuits, the rating is the sum of individual ratings divided by the number of circuits.
This provides a conservative estimate that accounts for potential overestimation of total capacity.
"""
function get_equivalent_rating(bp::BranchesParallel)
    # Sum of ratings divided by number of circuits
    return sum(PSY.get_rating(branch) for branch in bp.branches) / length(bp.branches)
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
    if isabstracttype(T)
        @warn "Parallel circuit contains mixed branch types, filters might be applied to more components than intended. Use Logging.Debug for additional information."
        @debug "Parallel circuit branch types: $(keys(double_circuit.branches))"
        @debug "Parallel circuit branch names: $(vcat([PSY.get_name.(v) for (k , v) in double_circuit.branches]))"
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
