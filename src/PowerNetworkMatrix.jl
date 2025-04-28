"""
The container code for PowerNetworkMatrix is based in JuMP's Container in order to
remove the limitations of AxisArrays and the doubts about long term maintenance
https://github.com/JuliaOpt/JuMP.jl/blob/master/src/Containers/DenseAxisArray.jl
JuMP'sCopyright 2017, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0.
"""

""" Type PowerNetworkMatrix gathers all the different types of Matrices considered in this package """
abstract type PowerNetworkMatrix{T} <: AbstractArray{T, 2} end

"""
Evaluates the bus indices for the given branch.

# Arguments
- `branch::PSY.ACTransmission`:
        system's branch
- `bus_lookup::Dict{Int, Int}`: A dictionary mapping bus numbers to their indices.
- `reverse_bus_search_map::Dict{Int, Int}`: map for reduced buses.

# Returns
- `fr_b::Int`:
        from bus index
- `to_b::Int`:
        to bus index
"""
function get_bus_indices(
    branch::PSY.ACTransmission,
    bus_lookup::Dict{Int, Int},
    reverse_bus_search_map::Dict{Int, Int},
)
    fr_bus_number = PSY.get_number(PSY.get_from(PSY.get_arc(branch)))
    to_bus_number = PSY.get_number(PSY.get_to(PSY.get_arc(branch)))
    if haskey(reverse_bus_search_map, fr_bus_number)
        fr_bus_ix = bus_lookup[reverse_bus_search_map[fr_bus_number]]
    else
        fr_bus_ix = bus_lookup[fr_bus_number]
    end
    if haskey(reverse_bus_search_map, to_bus_number)
        to_bus_ix = bus_lookup[reverse_bus_search_map[to_bus_number]]
    else
        to_bus_ix = bus_lookup[to_bus_number]
    end
    return fr_bus_ix, to_bus_ix
end

"""
    get_bus_indices(branch::PSY.Transformer3W, bus_lookup::Dict{Int, Int})

Retrieve the bus indices for the primary-secondary, secondary-tertiary, and primary-tertiary arcs of a 3-winding transformer.

# Arguments
- `branch::PSY.Transformer3W`: The 3-winding transformer branch.
- `bus_lookup::Dict{Int, Int}`: A dictionary mapping bus numbers to their indices.
- `reverse_bus_search_map::Dict{Int, Int}`: map for reduced buses.

# Returns
- A tuple containing three tuples:
  - `(ps_from, ps_to)`: from and to buses for the primary-secondary arc.
  - `(st_from, st_to)`: from and to buses for the secondary-tertiary arc.
  - `(pt_from, pt_to)`: from and to buses for the primary-tertiary arc.
"""
function get_bus_indices(
    branch::PSY.Transformer3W,
    bus_lookup::Dict{Int, Int},
    reverse_bus_search_map::Dict{Int, Int},
)
    ps_from_bus_number = PSY.get_number(PSY.get_from(PSY.get_primary_secondary_arc(branch)))
    if haskey(reverse_bus_search_map, ps_from_bus_number)
        ps_from = bus_lookup[reverse_bus_search_map[ps_from_bus_number]]
    else
        ps_from = bus_lookup[ps_from_bus_number]
    end

    ps_to_bus_number = PSY.get_number(PSY.get_to(PSY.get_primary_secondary_arc(branch)))
    if haskey(reverse_bus_search_map, ps_to_bus_number)
        ps_to = bus_lookup[reverse_bus_search_map[ps_to_bus_number]]
    else
        ps_to = bus_lookup[ps_to_bus_number]
    end

    st_from_bus_number =
        PSY.get_number(PSY.get_from(PSY.get_secondary_tertiary_arc(branch)))
    if haskey(reverse_bus_search_map, st_from_bus_number)
        st_from = bus_lookup[reverse_bus_search_map[st_from_bus_number]]
    else
        st_from = bus_lookup[st_from_bus_number]
    end

    st_to_bus_number = PSY.get_number(PSY.get_to(PSY.get_secondary_tertiary_arc(branch)))
    if haskey(reverse_bus_search_map, st_to_bus_number)
        st_to = bus_lookup[reverse_bus_search_map[st_to_bus_number]]
    else
        st_to = bus_lookup[st_to_bus_number]
    end

    pt_from_bus_number = PSY.get_number(PSY.get_from(PSY.get_primary_tertiary_arc(branch)))
    if haskey(reverse_bus_search_map, pt_from_bus_number)
        pt_from = bus_lookup[reverse_bus_search_map[pt_from_bus_number]]
    else
        pt_from = bus_lookup[pt_from_bus_number]
    end

    pt_to_bus_number = PSY.get_number(PSY.get_to(PSY.get_primary_tertiary_arc(branch)))
    if haskey(reverse_bus_search_map, pt_to_bus_number)
        pt_to = bus_lookup[reverse_bus_search_map[pt_to_bus_number]]
    else
        pt_to = bus_lookup[pt_to_bus_number]
    end
    return (ps_from, ps_to), (st_from, st_to), (pt_from, pt_to)
end

"""
Evaluates the map linking the system's buses and branches.

# Arguments
- `buses::AbstractVector{PSY.ACBus}`:
        system's buses
"""
function make_ax_ref(buses::AbstractVector{PSY.ACBus})
    return make_ax_ref(PSY.get_number.(buses))
end

"""
Checks if repetitions are present in the dictionary mapping buses and branches.

# Arguments
- `ax::AbstractVector`:
        generic abstract vector
"""
function make_ax_ref(ax::AbstractVector)
    ref = Dict{eltype(ax), Int}()
    for (ix, el) in enumerate(ax)
        if haskey(ref, el)
            @error("Repeated index element $el. Index sets must have unique elements.")
        end
        ref[el] = ix
    end
    return ref
end

# AbstractArray interface: overloadig methods
Base.isempty(A::PowerNetworkMatrix) = isempty(A.data)
Base.size(A::PowerNetworkMatrix) = size(A.data)
Base.LinearIndices(A::PowerNetworkMatrix) =
    error("PowerSystems.$(typeof(A)) does not support this operation.")
Base.axes(A::PowerNetworkMatrix) = A.axes
Base.CartesianIndices(A::PowerNetworkMatrix) =
    error("PowerSystems.$(typeof(A)) does not support this operation.")

# Indexing ###################################################################

Base.eachindex(A::PowerNetworkMatrix) = CartesianIndices(size(A.data))

"""
Gets bus indices to a certain branch index
"""
function lookup_index(i, lookup::Dict)
    return isa(i, Colon) ? Colon() : lookup[i]
end

"""
Gets bus indices to a certain branch name

# Arguments
- `i::PSY.ACBranch`:
        Power System AC branch
- `lookup::Dict`:
        Dictionary mapping branch and buses
"""
function lookup_index(i::PSY.ACBranch, lookup::Dict)
    return isa(i, Colon) ? Colon() : lookup[Base.to_index(i)]
end

"""
Gets bus indices to a certain branch name

# Arguments
- `i::PSY.ACBranch`:
        Power System AC branch
- `lookup::Dict`:
        Dictionary mapping branches and buses
"""
function lookup_index(i::PSY.ACBus, lookup::Dict)
    return isa(i, Colon) ? Colon() : lookup[Base.to_index(i)]
end

# Lisp-y tuple recursion trick to handle indexing in a nice type-
# stable way. The idea here is that `_to_index_tuple(idx, lookup)`
# performs a lookup on the first element of `idx` and `lookup`,
# then recurses using the remaining elements of both tuples.
# The compiler knows the lengths and types of each tuple, so
# all of the types are inferable.
function _to_index_tuple(idx::Tuple, lookup::Tuple)
    return tuple(
        lookup_index(first(idx), first(lookup)),
        _to_index_tuple(Base.tail(idx), Base.tail(lookup))...,
    )
end

# Handle the base case when we have more indices than lookups:
function _to_index_tuple(idx::NTuple{N}, ::NTuple{0}) where {N}
    return ntuple(k -> begin
            i = idx[k]
            (i == 1) ? 1 : error("invalid index $i")
        end, Val(N))
end

# Handle the base case when we have fewer indices than lookups:
_to_index_tuple(idx::NTuple{0}, lookup::Tuple) = ()

# Resolve ambiguity with the above two base cases
_to_index_tuple(idx::NTuple{0}, lookup::NTuple{0}) = ()

"""
Given the indices, gets the values of the power network matrix considered
"""
to_index(A::PowerNetworkMatrix, idx...) = _to_index_tuple(idx, A.lookup)

# Doing `Colon() in idx` is relatively slow because it involves
# a non-unrolled loop through the `idx` tuple which may be of
# varying element type. Another lisp-y recursion trick fixes that
has_colon(idx::Tuple{}) = false
has_colon(idx::Tuple) = isa(first(idx), Colon) || has_colon(Base.tail(idx))

# TODO: better error (or just handle correctly) when user tries to index with a range like a:b
# overloading other methods to consider PowerNetworkMatrix
function Base.getindex(A::PowerNetworkMatrix, row, column)
    i, j = to_index(A, row, column)
    return A.data[i, j]
end
Base.getindex(A::PowerNetworkMatrix, idx::CartesianIndex) = A.data[idx]
Base.setindex!(A::PowerNetworkMatrix, v, idx...) = A.data[to_index(A, idx...)...] = v
Base.setindex!(A::PowerNetworkMatrix, v, idx::CartesianIndex) = A.data[idx] = v

Base.IndexStyle(::Type{PowerNetworkMatrix}) = IndexAnyCartesian()

# Keys #######################################################################

"""
Structure to store the keys of a power network matrix

# Arguments
- `I<:Tuple`:
        turple containing the indices of the matrix
"""
struct PowerNetworkMatrixKey{T <: Tuple}
    I::T
end

Base.getindex(k::PowerNetworkMatrixKey, args...) = getindex(k.I, args...)

"""
Structure to store the keys of a power network matrix

# Arguments
- `product_iter::Base.Iterators.ProductIterator{T} where T <: Tuple`:
        iterator of the indices of the network power matrix
"""
struct PowerNetworkMatrixKeys{T <: Tuple}
    product_iter::Base.Iterators.ProductIterator{T}
end

# overloading methods
Base.length(iter::PowerNetworkMatrixKeys) = length(iter.product_iter)
function Base.eltype(iter::PowerNetworkMatrixKeys)
    return PowerNetworkMatrixKey{eltype(iter.product_iter)}
end
function Base.iterate(iter::PowerNetworkMatrixKeys)
    next = iterate(iter.product_iter)
    return next === nothing ? nothing : (PowerNetworkMatrixKey(next[1]), next[2])
end
function Base.iterate(iter::PowerNetworkMatrixKeys, state)
    next = iterate(iter.product_iter, state)
    return next === nothing ? nothing : (PowerNetworkMatrixKey(next[1]), next[2])
end
function Base.keys(a::PowerNetworkMatrix)
    return PowerNetworkMatrixKeys(Base.Iterators.product(a.axes...))
end
Base.getindex(a::PowerNetworkMatrix, k::PowerNetworkMatrixKey) = a[k.I...]

########
# Show #
########

# Adapted printing from JuMP's implementation of the Julia's show.jl
# used in PowerNetworkMatrixs

# Copyright (c) 2009-2016: Jeff Bezanson, Stefan Karpinski, Viral B. Shah,
# and other contributors:
#
# https://github.com/JuliaLang/julia/contributors
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.

function Base.summary(io::IO, A::PowerNetworkMatrix)
    _summary(io, A)
    for (k, ax) in enumerate(A.axes)
        print(io, "    Dimension $k, ")
        show(IOContext(io, :limit => true), ax)
        println(io)
    end
    print(io, "And data with size ", size(A))
    return
end

_summary(io::IO, ::T) where {T <: PowerNetworkMatrix} = println(io, "$T")

function Base.summary(
    io::IOContext{Base.GenericIOBuffer{Array{UInt8, 1}}},
    ::PowerNetworkMatrix,
)
    println(io, "PowerNetworkMatrix")
    return
end

function Base.summary(A::PowerNetworkMatrix)
    io = IOBuffer()
    Base.summary(io, A)
    String(take!(io))
    return
end

if isdefined(Base, :print_array) # 0.7 and later
    Base.print_array(io::IO, X::PowerNetworkMatrix) = Base.print_matrix(io, X.data)
end

# n-dimensional arrays
function Base.show_nd(
    io::IO,
    a::PowerNetworkMatrix,
    print_matrix::Function,
    label_slices::Bool,
)
    limit::Bool = get(io, :limit, false)
    if isempty(a)
        return
    end
    tailinds = Base.tail(Base.tail(axes(a.data)))
    nd = ndims(a) - 2
    for I in CartesianIndices(tailinds)
        idxs = I.I
        if limit
            for i in 1:nd
                ii = idxs[i]
                ind = tailinds[i]
                if length(ind) > 10
                    if ii == ind[4] && all(d -> idxs[d] == first(tailinds[d]), 1:(i - 1))
                        for j in (i + 1):nd
                            szj = size(a.data, j + 2)
                            indj = tailinds[j]
                            if szj > 10 && first(indj) + 2 < idxs[j] <= last(indj) - 3
                                @goto skip
                            end
                        end
                        #println(io, idxs)
                        print(io, "...\n\n")
                        @goto skip
                    end
                    if ind[3] < ii <= ind[end - 3]
                        @goto skip
                    end
                end
            end
        end
        if label_slices
            print(io, "[:, :, ")
            for i in 1:(nd - 1)
                show(io, a.axes[i + 2][idxs[i]])
                print(io, ", ")
            end
            show(io, a.axes[end][idxs[end]])
            println(io, "] =")
        end
        slice = view(a.data, axes(a.data, 1), axes(a.data, 2), idxs...)
        Base.print_matrix(io, slice)
        print(io, idxs == map(last, tailinds) ? "" : "\n\n")
        @label skip
    end
end

function Base.show(io::IO, array::PowerNetworkMatrix)
    summary(io, array)
    isempty(array) && return
    println(io, ":")
    Base.print_array(io, array)
    return
end

Base.to_index(b::PSY.ACBus) = PSY.get_number(b)
Base.to_index(b::T) where {T <: PSY.ACBranch} = PSY.get_name(b)

"""returns the raw array data of the `PowerNetworkMatrix`"""
get_data(mat::PowerNetworkMatrix) = mat.data

"""
    returns the lookup tuple of the `PowerNetworkMatrix`. The first entry corresponds
    to the first dimension and the second entry corresponds to the second dimension. For
    instance in Ybus the first dimension is buses and second dimension is buses too, and in
    PTDF the first dimension is branches and the second dimension is buses
"""
get_lookup(mat::PowerNetworkMatrix) = mat.lookup
