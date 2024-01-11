"""
Checks the network connectivity of the system using Depth First Search (DFS)
"""
function validate_connectivity(sys::PSY.System)
    sbn, _ = _find_subnetworks(sys)
    return length(sbn) == 1
end

function _find_subnetworks(sys::PSY.System)
    buses = get_buses(sys)
    branches = get_ac_branches(sys)
    @info "Validating connectivity with depth first search (network traversal)"
    M, bus_lookup = calculate_adjacency(branches, buses)
    slack_positions = find_slack_positions(buses, bus_lookup)
    return find_subnetworks(M, PSY.get_number.(buses)), slack_positions
end

"""
Finds the subnetworks in a system using Depth First Search (DFS). Returns a dictionary keyed
by the reference bus of the subnetworks if they exist
"""
function find_subnetworks(sys::PSY.System)
    sbn, ref_bus_positions = _find_subnetworks(sys)
    return assign_reference_buses(sbn, ref_bus_positions)
end
