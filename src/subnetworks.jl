"""
Takes the reference bus numbers and re-assigns the keys in the subnetwork dictionaries to use
the reference bus withing each subnetwork.
"""
function assign_reference_buses!(
    subnetworks::Dict{Int, Set{Int}},
    ref_buses::Set{Int},
)
    if isempty(ref_buses) || length(ref_buses) != length(subnetworks)
        @warn "The reference bus positions are not consistent with the subnetworks. References buses will be assigned arbitrarily"
        return deepcopy(subnetworks)
    end
    bus_groups = Dict{Int, Set{Int}}()
    for (bus_key, subnetwork_buses) in subnetworks
        ref_bus = intersect(ref_buses, subnetwork_buses)
        if length(ref_bus) == 1
            bus_groups[first(ref_bus)] = pop!(subnetworks, bus_key)
        elseif length(ref_bus) == 0
            @warn "No reference bus in the subnetwork associated with bus $bus_key. Reference bus assigned arbitrarily"
        elseif length(ref_bus) > 1
            error(
                "More than one reference bus in the subnetwork associated with bus $bus_key",
            )
        else
            @assert false
        end
    end
    return bus_groups
end

function assign_reference_buses!(
    subnetworks::Dict{Int, Set{Int}},
    ref_bus_positions::Set{Int},
    bus_lookup::Dict{Int, Int},
)
    ref_buses = [k for (k, v) in bus_lookup if v in ref_bus_positions]
    return assign_reference_buses!(subnetworks, ref_buses)
end
