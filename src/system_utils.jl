"""
Checks the network connectivity of the system.

# Keyword arguments
- `connectivity_method::Function = goderya_connectivity`: Specifies the method used as Goderya's algorithm (`goderya_connectivity`) or depth first search/network traversal (`dfs_connectivity`)
* Note that the default Goderya method is more efficient, but is resource intensive and may not scale well on large networks.
"""
function validate_connectivity(sys::PSY.System)
    sbn, _ = _find_subnetworks(sys)
    return length(sbn) == 1
end

function _find_subnetworks(sys::PSY.System)
    nodes = sort!(
        collect(
            PSY.get_components(x -> PSY.get_bustype(x) != BusTypes.ISOLATED, PSY.Bus, sys),
        );
        by = x -> PSY.get_number(x),
    )
    branches = get_ac_branches(sys)
    @info "Validating connectivity with depth first search (network traversal)"
    M, bus_lookup = calculate_adjacency(branches, nodes)
    return find_subnetworks(M, collect(keys(bus_lookup))), find_slack_positions(nodes)
end

function find_subnetworks(sys::PSY.System)
    sbn, ref_bus_positions = _find_subnetworks(sys)
    return assing_reference_buses(sbn, ref_bus_positions)
end
