"""
Incidence matrix: shows connection between buses, defining arcs

# Arguments
- `data::SparseArrays.SparseMatrixCSC{Int8, Int}`:
        the actual Incidence matrix.
- `axes<:NTuple{2, Dict}`:
        Tuple containing two vectors (the first one showing the branches names,
        the second showing the buses numbers).
- `lookup<:NTuple{2, Dict}`:
        Tuple containing two dictionaries, the first mapping the branches
        and buses with their enumerated indexes.
- `ref_bus_positions::Set{Int}`:
        Vector containing the indices of the reference slack buses.
- `network_reduction::NetworkReduction`:
        Structure containing the details of the network reduction applied when computing the matrix
"""
struct IncidenceMatrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Int8}
    data::SparseArrays.SparseMatrixCSC{Int8, Int}
    axes::Ax
    lookup::L
    ref_bus_numbers::Set{Int}
    network_reduction::NetworkReduction
end

# functions to get stored data
get_axes(A::IncidenceMatrix) = A.axes
get_lookup(A::IncidenceMatrix) = A.lookup
get_slack_positions(A::IncidenceMatrix) = [A.lookup[x] for x in A.ref_bus_numbers]

function IncidenceMatrix(sys::PSY.System;
    check_connectivity::Bool = true,
    network_reductions::Vector{NetworkReductionTypes} = NetworkReductionTypes[],
    kwargs...,
)
    return IncidenceMatrix(
        Ybus(
            sys;
            check_connectivity = check_connectivity,
            network_reductions = network_reductions,
            kwargs...,
        ),
    )
end

function IncidenceMatrix(ybus::Ybus)
    nr = ybus.network_reduction
    direct_arcs = [x for x in keys(nr.direct_branch_map)]
    parallel_arcs = [x for x in keys(nr.parallel_branch_map)]
    series_arcs = [x for x in keys(nr.series_branch_map)]
    transformer_arcs = [x for x in keys(nr.transformer3W_map)]
    bus_ax = ybus.axes[1]
    bus_lookup = ybus.lookup[1]
    arc_ax = vcat(direct_arcs, parallel_arcs, series_arcs, transformer_arcs)
    n_entries = length(arc_ax) * 2
    A_I = Vector{Int}(undef, n_entries)
    A_J = Vector{Int}(undef, n_entries)
    A_V = Vector{Int8}(undef, n_entries)
    for (ix, arc) in enumerate(arc_ax)
        A_I[2 * ix - 1] = ix
        A_J[2 * ix - 1] = get_bus_index(arc[1], bus_lookup, nr)
        A_V[2 * ix - 1] = 1
        A_I[2 * ix] = ix
        A_J[2 * ix] = get_bus_index(arc[2], bus_lookup, nr)
        A_V[2 * ix] = -1
    end
    data = SparseArrays.sparse(A_I, A_J, A_V)
    axes = (arc_ax, bus_ax)
    lookup = (make_ax_ref(arc_ax), make_ax_ref(bus_ax))
    ref_bus_numbers = ybus.ref_bus_numbers
    return IncidenceMatrix(data, axes, lookup, ref_bus_numbers, ybus.network_reduction)
end
