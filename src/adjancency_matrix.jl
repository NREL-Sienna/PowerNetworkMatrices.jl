
"""
Nodal incidence matrix (Adjacency) is an N x N matrix describing a power system with N buses. It represents the directed connectivity of the buses in a power system.

The AdjacencyMatrix Struct is indexed using the Bus Numbers, no need for them to be sequential
"""
struct AdjacencyMatrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Int8}
    data::SparseArrays.SparseMatrixCSC{Int8, Int}
    axes::Ax
    lookup::L
    ref_bus_positions::Vector{Int}
end

# functions to get stored data
get_axes(A::AdjacencyMatrix) = A.axes
get_lookup(A::AdjacencyMatrix) = A.lookup
get_slack_position(A::AdjacencyMatrix) = A.ref_bus_positions

"""
Builds a AdjacencyMatrix from the system. The return is an N x N AdjacencyMatrix Array indexed with the bus numbers.

# Keyword arguments
- `check_connectivity::Bool`: Checks connectivity of the network using Goderya's algorithm
- `connectivity_method::Function = goderya_connectivity`: method (`goderya_connectivity` or `dfs_connectivity`) for connectivity validation
"""
function AdjacencyMatrix(sys::PSY.System; check_connectivity::Bool = true, kwargs...)
    nodes = sort!(
        collect(
            PSY.get_components(x -> PSY.get_bustype(x) != BusTypes.ISOLATED, PSY.Bus, sys),
        );
        by = x -> PSY.get_number(x),
    )
    branches = PSY.get_components(PSY.get_available, PSY.Branch, sys)
    return AdjacencyMatrix(
        branches,
        nodes;
        check_connectivity = check_connectivity,
        kwargs...,
    )
end

"""
Builds a AdjacencyMatrix from a collection of buses and branches. The return is an N x N AdjacencyMatrix Array indexed with the bus numbers.

# Keyword arguments
- `check_connectivity::Bool`: Checks connectivity of the network using Goderya's algorithm
"""
function AdjacencyMatrix(
    branches,
    nodes::Vector{PSY.Bus};
    check_connectivity::Bool = true,
    kwargs...,
)
    a, bus_lookup = calculate_adjacency(branches, nodes)
    bus_ax = PSY.get_number.(nodes)
    axes = (bus_ax, bus_ax)
    look_up = (bus_lookup, bus_lookup)

    if check_connectivity
        connected = _validate_connectivity(a, bus_lookup, dfs_connectivity)
        !connected && throw(DataFormatError("Network not connected"))
    end

    return AdjacencyMatrix(a, axes, look_up, find_slack_positions(nodes))
end

function validate_connectivity(
    M::AdjacencyMatrix;
    connectivity_method::Function = goderya_connectivity,
)
    return _validate_connectivity(M.data, M.lookup[1], connectivity_method)
end

function _validate_connectivity(
    M::SparseArrays.SparseMatrixCSC{Int8, Int},
    bus_lookup::Dict{Int, Int},
    connectivity_method::Function,
)
    connected = connectivity_method(M, bus_lookup)
    return connected
end

function _goderya(ybus::SparseArrays.SparseMatrixCSC)
    node_count = size(ybus)[1]
    max_I = node_count^2
    I, J, val = SparseArrays.findnz(ybus)
    T = SparseArrays.sparse(I, J, ones(Int8, length(val)))
    T_ = T * T
    for n in 1:(node_count - 1)
        I, _, _ = SparseArrays.findnz(T_)
        if length(I) == max_I
            @info "The System has no islands"
            break
        elseif length(I) < max_I
            temp = T_ * T
            I_temp, _, _ = SparseArrays.findnz(temp)
            if all(I_temp == I)
                @warn "The system contains islands" maxlog = 1
            end
            T_ = temp
        else
            error("Invalid connectivity condition")
        end
        #@assert n < node_count - 1
    end
    return I
end

function goderya_connectivity(M::SparseArrays.SparseMatrixCSC{Int8, Int},
    bus_lookup::Dict{Int, Int})
    @info "Validating connectivity with Goderya algorithm"
    node_count = length(bus_lookup)

    if node_count > GODERYA_MAX_PERFORMANCE_NODE
        @warn "The Goderya algorithm is memory intensive on large networks and may not scale well, try `connectivity_method = dfs_connectivity"
    end

    I = _goderya(M)

    connections = Dict(i => count(x -> x == i, I) for i in Set(I))

    if length(Set(I)) == node_count
        connected = true
        if any(values(connections) .!= node_count)
            cc = Set(values(connections))
            @warn "Network has at least $(length(cc)) connected components with $cc nodes"
            connected = false
        end
    else
        connected = false
    end
    return connected
end

function dfs_connectivity(M::SparseArrays.SparseMatrixCSC,
    bus_lookup::Dict{Int, Int})
    @info "Validating connectivity with depth first search (network traversal)"
    cc = find_connected_components(M, bus_lookup)
    if length(cc) != 1
        @warn "Network has at least $(length(cc)) connected components with $(length.(cc)) nodes"
        connected = false
    else
        connected = true
    end
    return connected
end

function find_subnetworks(M::AdjacencyMatrix)
    bus_numbers = M.axes[2]
    return find_subnetworks(M.data, bus_numbers)
end
