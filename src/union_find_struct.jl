struct UnionFind
    rep::Dict{Int, Int} # maps each element to its representative
end

UnionFind(elems::Vector{Int}) = UnionFind(Dict{Int, Int}(x => x for x in elems))

function get_rep(uf::UnionFind, x::Int)
    while uf.rep[x] != x
        uf.rep[x] = uf.rep[uf.rep[x]] # path compression
        x = uf.rep[x]
    end
    return x
end

function union_sets!(uf::UnionFind, x::Int, y::Int)
    rootX = get_rep(uf, x)
    rootY = get_rep(uf, y)
    if rootX != rootY
        uf.rep[rootY] = min(rootX, rootY)
        uf.rep[rootX] = min(rootX, rootY)
    end
end
