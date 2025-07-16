"""
Checks the network connectivity of the system using Depth First Search (DFS)
"""
function validate_connectivity(sys::PSY.System)
    sbn, _ = _find_subnetworks(sys)
    return length(sbn) == 1
end

function _find_subnetworks(sys::PSY.System)
    Ymatrix = Ybus(sys; check_connectivity = false)
    @info "Validating connectivity with depth first search (network traversal)"
    subnetworks = find_subnetworks(Ymatrix)
    return subnetworks, Ymatrix.ref_bus_numbers
end

"""
Finds the subnetworks in a system using Depth First Search (DFS). Returns a dictionary keyed
by the reference bus of the subnetworks if they exist
"""
function find_subnetworks(sys::PSY.System)
    sbn, ref_buses = _find_subnetworks(sys)
    return assign_reference_buses!(sbn, ref_buses)
end
