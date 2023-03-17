"""
Finds the set of bus numbers that belong to each connnected component in the System
"""
function find_connected_components(sys::PSY.System)
    a = Adjacency(sys; check_connectivity = false)
    return find_connected_components(a.data, a.lookup[1])
end

"""
Checks the network connectivity of the system.

# Keyword arguments
- `connectivity_method::Function = goderya_connectivity`: Specifies the method used as Goderya's algorithm (`goderya_connectivity`) or depth first search/network traversal (`dfs_connectivity`)
* Note that the default Goderya method is more efficient, but is resource intensive and may not scale well on large networks.
"""
function validate_connectivity(
    sys::PSY.System;
    connectivity_method::Function = goderya_connectivity,
)
    nodes = sort!(
        collect(
            PSY.get_components(x -> PSY.get_bustype(x) != BusTypes.ISOLATED, PSY.Bus, sys),
        );
        by = x -> PSY.get_number(x),
    )
    branches = get_ac_branches(sys)
    a = AdjacencyMatrix(branches, nodes; check_connectivity = false)

    return validate_connectivity(
        a;
        connectivity_method = connectivity_method,
    )
end
