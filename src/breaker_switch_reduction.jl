"""
Builds a NetworkReduction by removing branches with zero impedance. 

# Arguments
- `sys::System`
"""
function get_breaker_switch_reduction(sys::PSY.System; prior_reduction = NetworkReduction())
    validate_reduction_type(
        NetworkReductionTypes.BREAKER_SWITCH,
        get_reduction_type(prior_reduction),
    )
    retained_branches = Set{String}()
    removed_branches = Set{String}()
    A = IncidenceMatrix(sys)
    bus_map = A.lookup[2]
    bus_reduction_map_index = Dict{Int, Set{Int}}(k => Set{Int}() for k in keys(bus_map))
    for br in get_ac_branches(sys, prior_reduction.removed_branches)
        if isa(br, PSY.DiscreteControlledACBranch)
            r = PSY.get_r(br)
            x = PSY.get_x(br)
            status = PSY.get_branch_status(br)
            if r == 0.0 && x < ZERO_IMPEDANCE_LINE_REACTANCE_THRESHOLD
                if status == PSY.DiscreteControlledBranchStatus.CLOSED
                    from_bus_number = PSY.get_number(PSY.get_from(PSY.get_arc(br)))
                    to_bus_number = PSY.get_number(PSY.get_to(PSY.get_arc(br)))
                    from_bus_reduced = !haskey(bus_reduction_map_index, from_bus_number)
                    to_bus_reduced = !haskey(bus_reduction_map_index, to_bus_number)
                    if from_bus_reduced && !to_bus_reduced
                        pop!(bus_reduction_map_index, to_bus_number)
                        reduction_set = Set{Int}(to_bus_number)
                        new_number = _find_reduced(bus_reduction_map_index, from_bus_number)
                        union!(bus_reduction_map_index[new_number], reduction_set)
                    elseif !from_bus_reduced && to_bus_reduced
                        pop!(bus_reduction_map_index, from_bus_number)
                        reduction_set = Set{Int}(from_bus_number)
                        new_number = _find_reduced(bus_reduction_map_index, to_bus_number)
                        union!(bus_reduction_map_index[new_number], reduction_set)
                    elseif !from_bus_reduced && !to_bus_reduced
                        pop!(bus_reduction_map_index, from_bus_number)
                        reduction_set = Set{Int}(from_bus_number)
                        union!(bus_reduction_map_index[to_bus_number], reduction_set)
                    end
                end
                push!(removed_branches, PSY.get_name(br))
            else
                push!(retained_branches, PSY.get_name(br))
            end
        end
    end
    reverse_bus_search_map =
        _make_reverse_bus_search_map(bus_reduction_map_index, length(bus_map))
    new_reduction = NetworkReduction(
        bus_reduction_map_index,
        reverse_bus_search_map,
        removed_branches,
        retained_branches,
        Vector{PSY.ACTransmission}(),
        Vector{PSY.FixedAdmittance}(),
        [NetworkReductionTypes.BREAKER_SWITCH],
    )
    if isempty(prior_reduction)
        return new_reduction
    else
        return compose_reductions(prior_reduction, new_reduction, length(bus_map))
    end
end

function _find_reduced(bus_reduction_map_index, number)
    for (k, v) in bus_reduction_map_index
        if number ∈ v
            return k
        end
    end
end
