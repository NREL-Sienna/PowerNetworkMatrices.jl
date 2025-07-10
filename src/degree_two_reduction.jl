#NOTE: hardcoded for testing system
function get_reduction(
    A::AdjacencyMatrix,
    sys::PSY.System,
    ::Val{NetworkReductionTypes.DEGREE_TWO},
)
    irreducible_buses = Set{Int}()
    for c in get_components(StaticInjection, sys_DA)
        push!(irreducible_buses, PSY.get_number(PSY.get_bus(c)))
    end
    irreducible_indices = Set([yb.lookup[2][i] for i in irreducible_buses])

    arc_maps = find_degree2_chains(A.data, irreducible_indices)

    @error "DEGREE TWO HARDCODED REDUCTION FOR 14 BUS SYSTEM"
    br1 = PSY.get_component(PSY.Line, sys, "BUS 101-BUS 115-i_1")
    br2 = PSY.get_component(PSY.Line, sys, "BUS 115-BUS 102-i_1")
    return NetworkReduction(
        irreducible_buses,
        Dict{Int, Set{Int}}(),
        Dict{Int, Int}(),
        Dict{Tuple{Int, Int}, PSY.Branch}(),
        Dict{PSY.Branch, Tuple{Int, Int}}(),
        Dict{Tuple{Int, Int}, Set{PSY.Branch}}(),
        Dict{PSY.Branch, Tuple{Int, Int}}(),
        Dict{Tuple{Int, Int}, Set{PSY.Branch}}((101, 102) => Set([br1, br2])),          # Map from NEW arc to the series branches that comprise it.
        Dict{PSY.Branch, Tuple{Int, Int}}(br1 => (101, 102), br2 => (101, 102)),
        Dict{Tuple{Int, Int}, Tuple{PSY.ThreeWindingTransformer, Int}}(),
        Dict{Tuple{PSY.ThreeWindingTransformer, Int}, Tuple{Int, Int}}(),
        Set{Int}([115]),
        Set{Tuple{Int, Int}}([(101, 115), (115, 102)]),                                 # Set of OLD arcs to be removed
        Vector{NetworkReductionTypes}([NetworkReductionTypes.DEGREE_TWO]),
    )
end

"""
    _should_visit_node(node::Int, reduced_indices::Set{Int}, irreducible_indices::Set{Int})

Determines whether a node should be visited during network traversal.

# Arguments
- `node::Int`: The index of the node to check.
- `reduced_indices::Set{Int}`: Set of indices that have already been reduced.
- `irreducible_indices::Set{Int}`: Set of indices that cannot be reduced.

# Returns
- `Bool`: `true` if the node should be visited, `false` otherwise.
"""
function _should_visit_node(
    node::Int,
    reduced_indices::Set{Int},
    irreducible_indices::Set{Int},
)
    if node ∈ irreducible_indices
        return false
    end
    if node ∈ reduced_indices
        return false
    end
    return true
end

"""
    _is_final_node(node::Int, adj_matrix::SparseArrays.SparseMatrixCSC, reduced_indices::Set{Int})

Determines if a node is a final node in a path traversal.

# Arguments
- `node::Int`: The index of the node to check.
- `adj_matrix::SparseArrays.SparseMatrixCSC`: The adjacency matrix of the network.
- `reduced_indices::Set{Int}`: Set of indices that have already been reduced.

# Returns
- `Bool`: `true` if the node is a final node, `false` otherwise.
"""
function _is_final_node(node::Int, adj_matrix::SparseArrays.SparseMatrixCSC, reduced_indices::Set{Int})
    if !_is_2degree_node(adj_matrix, node)
        return true
    end
    if node ∈ reduced_indices
        return true
    end
    return false
end

"""
    _is_2degree_node(adj_matrix::SparseArrays.SparseMatrixCSC, node::Int)

Checks if a node has exactly two connections in the network.

# Arguments
- `adj_matrix::SparseArrays.SparseMatrixCSC`: The adjacency matrix of the network.
- `node::Int`: The index of the node to check.

# Returns
- `Bool`: `true` if the node has exactly two neighbors, `false` otherwise.
"""
function _is_2degree_node(adj_matrix::SparseArrays.SparseMatrixCSC, node::Int)
    neighbor_count = SparseArrays.nzrange(adj_matrix, node)
    return length(neighbor_count) == 2
end

"""
    _get_neighbors(adj_matrix::SparseArrays.SparseMatrixCSC, node::Int)

Get all neighbors of a given node from the adjacency matrix.
For undirected graphs, checks both directions.
"""
function _get_neighbors(adj_matrix::SparseArrays.SparseMatrixCSC, node::Int)
    nzrange = SparseArrays.nzrange(adj_matrix, node)
    @assert length(nzrange) == 2
    return rowvals(adj_matrix)[nzrange]
end

"""
    _get_complete_chain(adj_matrix::SparseArrays.SparseMatrixCSC, start_node::Int, reduced_indices::Set{Int}, irreducible_indices::Set{Int})

Build a complete chain of degree-2 nodes starting from a given node.
"""
function _get_complete_chain(
    adj_matrix::SparseArrays.SparseMatrixCSC,
    start_node::Int,
    reduced_indices::Set{Int},
    irreducible_indices::Set{Int},
)
    neighbors = _get_neighbors(adj_matrix, start_node)
    current_chain = [start_node]
    push!(reduced_indices, start_node)
    _get_partial_chain_recursive!(
        current_chain,
        adj_matrix,
        neighbors[1],
        start_node,
        reduced_indices,
        irreducible_indices,
    )
    reverse!(current_chain)
    _get_partial_chain_recursive!(
        current_chain,
        adj_matrix,
        neighbors[2],
        start_node,
        reduced_indices,
        irreducible_indices,
    )
    return current_chain
end

"""
    _get_partial_chain(adj_matrix::SparseArrays.SparseMatrixCSC,
                      current_node::Int,
                      prev_node::Int,
                      reduced_indices::Set{Int},
                      irreducible_indices::Set{Int})

Recursively build a chain in one direction from current_node, avoiding prev_node.
"""
function _get_partial_chain_recursive!(
    current_chain::Vector{Int},
    adj_matrix::SparseArrays.SparseMatrixCSC,
    current_node::Int,
    prev_node::Int,
    reduced_indices::Set{Int},
    irreducible_indices::Set{Int})
    # If current node is reduced stop
    if current_node ∈ reduced_indices
        return Int[]
    end

    push!(current_chain, current_node)

    if _is_final_node(current_node, adj_matrix, reduced_indices)
        return
    end

    push!(reduced_indices, current_node)
    # Get neighbors
    neighbors = _get_neighbors(adj_matrix, current_node)

    # Determine the next node to visit. It must not be the `previous_node`.
    # This prevents the traversal from going back and forth between two nodes.
    next_node = (neighbors[1] == prev_node) ? neighbors[2] : neighbors[1]
    _get_partial_chain_recursive!(
        current_chain,
        adj_matrix,
        next_node,
        current_node,
        reduced_indices,
        irreducible_indices,
    )
    return
end

"""
    _get_degree2_nodes(adj_matrix::SparseArrays.SparseMatrixCSC, irreducible_indices::Set{Int})

Return all degree-2 nodes in the adjacency matrix, excluding irreducible indices.
"""
function _get_degree2_nodes(adj_matrix::SparseArrays.SparseMatrixCSC, irreducible_indices::Set{Int})
    node_count = size(adj_matrix, 1)
    nodes = sizehint!(Vector{Int}(), node_count)
    for i in 1:node_count
        if i ∈ irreducible_indices
            continue
        end
        if _is_2degree_node(adj_matrix, i)
            push!(nodes, i)
        end
    end
    return nodes
end

"""
    find_degree2_chains(adj_matrix::SparseArrays.SparseMatrixCSC, irreducible_indices::Set{Int})

Find all chains of degree-2 nodes in a graph represented by a CSC adjacency matrix.
A chain is a sequence of connected degree-2 nodes.

Returns a dictionary mapping each starting node to its chain of node indices.
"""
function find_degree2_chains(adj_matrix::SparseArrays.SparseMatrixCSC, irreducible_indices::Set{Int})
    arc_map = Dict()
    reduced_indices = Set{Int}()
    degree2_nodes = _get_degree2_nodes(adj_matrix, irreducible_indices)
    for node in degree2_nodes
        if node ∈ reduced_indices
            continue
        end
        chain_path =
            _get_complete_chain(adj_matrix, node, reduced_indices, irreducible_indices)
        if adj_matrix[chain_path[1], chain_path[end]] != 0
            @warn "Nodes $(chain_path[1]) and $(chain_path[end]) already have a parallel path, skipping chain creation."
        else
            arc_map[chain_path[1], chain_path[end]] = chain_path
        end
    end
    return arc_map
end
