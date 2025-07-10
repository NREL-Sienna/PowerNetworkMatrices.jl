function _goderya(ybus::SparseArrays.SparseMatrixCSC)
    node_count = size(ybus)[1]
    max_I = node_count^2
    I, J, val = SparseArrays.findnz(ybus)
    T = SparseArrays.sparse(I, J, ones(Int, length(val)))
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
            @assert false
        end
        #@assert n < node_count - 1
    end
    return I
end

function validate_connectivity(
    M,
    nodes::Vector{PSY.ACBus},
    bus_lookup::Dict{Int64, Int64};
    connectivity_method::Function = goderya_connectivity,
)
    connected = connectivity_method(M, nodes, bus_lookup)
    return connected
end

function goderya_connectivity(M, nodes::Vector{PSY.ACBus}, bus_lookup::Dict{Int64, Int64})
    @info "Validating connectivity with Goderya algorithm"
    length(nodes) > 15_000 &&
        @warn "The Goderya algorithm is memory intensive on large networks and may not scale well, try `connectivity_method = dfs_connectivity"

    I = _goderya(M)

    node_count = length(nodes)
    connections = Dict([i => count(x -> x == i, I) for i in Set(I)])

    if length(Set(I)) == node_count
        connected = true
        if any(values(connections) .!= node_count)
            cc = Set(values(connections))
            @warn "Network has at least $(length(cc)) connected components with $cc nodes"
            connected = false
        end
    else
        disconnected_nodes = PSY.get_name.(nodes[setdiff(values(bus_lookup), I)])
        @warn "Principal connected component does not contain:" disconnected_nodes
        connected = false
    end
    return connected
end

"""
Finds the set of bus numbers that belong to each connected component in the System
"""
# this function extends the PowerModels.jl implementation to accept a System
function find_connected_components(sys::PSY.System)
    a = Adjacency(sys; check_connectivity = false)
    return find_connected_components(a.data, a.lookup[1])
end

# this function extends the PowerModels.jl implementation to accept an adjacency matrix and bus lookup
function find_connected_components(M, bus_lookup::Dict{Int64, Int64})
    pm_buses = Dict([i => Dict("bus_type" => 1, "bus_i" => b) for (i, b) in bus_lookup])

    arcs = findall((LinearAlgebra.UpperTriangular(M) - LinearAlgebra.I) .!= 0)
    pm_branches = Dict([
        i => Dict("f_bus" => a[1], "t_bus" => a[2], "br_status" => 1) for
        (i, a) in enumerate(arcs)
    ],)

    data = Dict("bus" => pm_buses, "branch" => pm_branches)
    cc = PSY.calc_connected_components(data)
    bus_decode = Dict(value => key for (key, value) in bus_lookup)
    connected_components = Vector{Set{Int64}}()
    for c in cc
        push!(connected_components, Set([bus_decode[b] for b in c]))
    end
    return Set(connected_components)
end

function dfs_connectivity(M, ::Vector{PSY.ACBus}, bus_lookup::Dict{Int64, Int64})
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

"""Find part of the union-find disjoint set data structure. Vector because nodes are 1:n."""
function get_representative(uf::Vector{Int}, x::Int)
    while uf[x] != x
        uf[x] = uf[uf[x]] # path compression
        x = uf[x]
    end
    return x
end

"""Union part of the union-find disjoint set data structure. Vector because nodes are 1:n."""
function union_sets!(uf::Vector{Int}, x::Int, y::Int)
    x == y && return
    rootX = get_representative(uf, x)
    rootY = get_representative(uf, y)
    if rootX != rootY
        uf[rootY] = rootX
    end
end

"""
Finds the subnetworks present in the considered System. This is evaluated by taking
a the ABA or Adjacency Matrix.

# Arguments
- `M::SparseArrays.SparseMatrixCSC`:
        input sparse matrix.
- `bus_numbers::Vector{Int}`:
        vector containing the indices of the system's buses.
"""
function find_subnetworks(M::SparseArrays.SparseMatrixCSC, bus_numbers::Vector{Int})
    rows = SparseArrays.rowvals(M)
    # find connected components working with indices, so can use a vector instead of a set.
    # uf for union-find data structure: initially, each bus as its own set of a single element.
    uf = collect(1:size(bus_numbers, 1))
    for ix in 1:size(bus_numbers, 1)
        neighbors = SparseArrays.nzrange(M, ix)
        if length(neighbors) <= 1
            @warn "Bus $(bus_numbers[ix]) is islanded"
            continue
        end
        for j in neighbors
            row_ix = rows[j]
            union_sets!(uf, ix, row_ix)
        end
    end
    for i in 1:length(uf)
        uf[i] = get_representative(uf, i)
    end
    # now we have the representatives, so assemble the subnetworks.
    subnetworks = Dict{Int, Set{Int}}()
    for (representative, bus_num) in zip(uf, bus_numbers)
        subnetwork = get!(subnetworks, bus_numbers[representative], Set{Int}())
        push!(subnetwork, bus_num)
    end
    return subnetworks
end
