@kwdef struct DegreeTwoReduction <: NetworkReduction
    irreducible_buses::Vector{Int} = Vector{Int}()
    reduce_reactive_power_injectors::Bool = true
end
get_irreducible_buses(nr::DegreeTwoReduction) = nr.irreducible_buses
get_reduce_reactive_power_injectors(nr::DegreeTwoReduction) =
    nr.reduce_reactive_power_injectors

function get_degree2_reduction(
    A::AdjacencyMatrix,
    sys::PSY.System,
    irreducible_buses::Set{Int},
    reduction::DegreeTwoReduction,
)
    network_reduction_data = A.network_reduction_data
    reverse_bus_search_map = network_reduction_data.reverse_bus_search_map
    for c in PSY.get_components(PSY.StaticInjection, sys)
        push!(irreducible_buses, PSY.get_number(PSY.get_bus(c)))
    end
    irreducible_non_reduced_buses =
        unique([get(reverse_bus_search_map, k, k) for k in irreducible_buses])
    irreducible_indices = Set([A.lookup[2][i] for i in irreducible_non_reduced_buses])
    reverse_bus_lookup = Dict(v => k for (k, v) in A.lookup[2])
    arc_maps = find_degree2_chains(A.data, irreducible_indices)
    series_branch_map = Dict{Tuple{Int, Int}, Vector{Any}}()
    reverse_series_branch_map = Dict{PSY.Branch, Tuple{Int, Int}}()
    removed_buses = Set{Int}()
    removed_arcs = Set{Tuple{Int, Int}}()
    for (composite_arc_ix, segment_ix) in arc_maps
        composite_arc = (
            reverse_bus_lookup[composite_arc_ix[1]],
            reverse_bus_lookup[composite_arc_ix[2]],
        )
        segment_numbers = [reverse_bus_lookup[x] for x in segment_ix]
        @assert composite_arc[1] == segment_numbers[1]
        @assert composite_arc[2] == segment_numbers[end]
        segments = Vector{Any}()
        for ix in 1:(length(segment_numbers) - 1)
            segment_arc = (segment_numbers[ix], segment_numbers[ix + 1])
            segment_arc, entry = _get_branch_map_entry(network_reduction_data, segment_arc)
            push!(segments, entry)
            push!(removed_arcs, segment_arc)
            ix != 1 && push!(removed_buses, segment_numbers[ix])
        end
        series_branch_map[composite_arc] = segments
    end

    #TODO - add the reverse_series_branch_map

    return NetworkReductionData(;
        irreducible_buses = irreducible_buses,
        series_branch_map = series_branch_map,
        reverse_series_branch_map = reverse_series_branch_map,
        removed_buses = removed_buses,
        removed_arcs = removed_arcs,
        reductions = NetworkReduction[reduction],
    )
end

function _get_branch_map_entry(nr::NetworkReductionData, arc::Tuple{Int, Int})
    reverse_arc = (arc[2], arc[1])
    direct_branch_map = nr.direct_branch_map
    parallel_branch_map = nr.parallel_branch_map
    series_branch_map = nr.series_branch_map
    transformer3W_map = nr.transformer3W_map

    if haskey(direct_branch_map, arc)
        return arc, direct_branch_map[arc]
    elseif haskey(direct_branch_map, reverse_arc)
        return reverse_arc, direct_branch_map[reverse_arc]
    elseif haskey(parallel_branch_map, arc)
        return arc, parallel_branch_map[arc]
    elseif haskey(parallel_branch_map, reverse_arc)
        return reverse_arc, parallel_branch_map[reverse_arc]
    elseif haskey(series_branch_map, arc)
        return arc, series_branch_map[arc]
    elseif haskey(series_branch_map, reverse_arc)
        return reverse_arc, series_branch_map[reverse_arc]
    elseif haskey(transformer3W_map, arc)
        return arc, transformer3W_map[arc]
    elseif haskey(transformer3W_map, reverse_arc)
        return reverse_arc, transformer3W_map[reverse_arc]
    else
        error("Arc $segment_arc not found in the existing network reduction mappings.")
    end
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
function _is_final_node(
    node::Int,
    adj_matrix::SparseArrays.SparseMatrixCSC,
    reduced_indices::Set{Int},
)
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
function _get_degree2_nodes(
    adj_matrix::SparseArrays.SparseMatrixCSC,
    irreducible_indices::Set{Int},
)
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
function find_degree2_chains(
    adj_matrix::SparseArrays.SparseMatrixCSC,
    irreducible_indices::Set{Int},
)
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
