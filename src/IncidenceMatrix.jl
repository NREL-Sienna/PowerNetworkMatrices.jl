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
    subnetwork_axes::Dict{Int, Ax}
    network_reduction_data::NetworkReductionData
end

# functions to get stored data
get_axes(M::IncidenceMatrix) = M.axes
get_lookup(M::IncidenceMatrix) = M.lookup
get_ref_bus(M::IncidenceMatrix) = collect(keys(M.subnetwork_axes))
get_ref_bus_position(M::IncidenceMatrix) =
    [get_bus_lookup(M)[x] for x in keys(M.subnetwork_axes)]
get_network_reduction_data(M::IncidenceMatrix) = M.network_reduction_data
get_arc_axis(M::IncidenceMatrix) = M.axes[1]
get_arc_lookup(M::IncidenceMatrix) = M.lookup[1]
get_bus_axis(M::IncidenceMatrix) = M.axes[2]
get_bus_lookup(M::IncidenceMatrix) = M.lookup[2]

function get_reduction(
    A::IncidenceMatrix,
    ::PSY.System,
    reduction::RadialReduction,
)
    irreducible_buses = get_irreducible_buses(reduction)
    irreducible_positions = Set([A.lookup[2][x] for x in irreducible_buses])
    exempt_positions = union(get_ref_bus_position(A), irreducible_positions)
    bus_reduction_map, reverse_bus_search_map, radial_arcs =
        calculate_radial_arcs(A.data, A.lookup[1], A.lookup[2], Set(exempt_positions))

    return NetworkReductionData(;
        irreducible_buses = Set(irreducible_buses),
        bus_reduction_map = bus_reduction_map,
        reverse_bus_search_map = reverse_bus_search_map,
        removed_arcs = radial_arcs,
        reductions = ReductionContainer(; radial_reduction = reduction),
    )
end

function IncidenceMatrix(sys::PSY.System;
    network_reductions::Vector{NetworkReduction} = NetworkReduction[],
    kwargs...,
)
    return IncidenceMatrix(
        Ybus(
            sys;
            network_reductions = network_reductions,
            kwargs...,
        ),
    )
end

function IncidenceMatrix(ybus::Ybus)
    nr = ybus.network_reduction_data
    direct_arcs = [x for x in keys(nr.direct_branch_map)]
    parallel_arcs = [x for x in keys(nr.parallel_branch_map)]
    series_arcs = [x for x in keys(nr.series_branch_map)]
    transformer_arcs = [x for x in keys(nr.transformer3W_map)]
    additional_arcs = [x for x in keys(nr.added_branch_map)]
    bus_ax = get_bus_axis(ybus)
    bus_lookup = ybus.lookup[1]
    arc_ax = unique(
        vcat(direct_arcs, parallel_arcs, series_arcs, transformer_arcs, additional_arcs),
    )
    n_isolated_buses = length(get_isolated_buses(ybus))
    n_entries = length(arc_ax) * 2 + n_isolated_buses
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
    A_I[(end - n_isolated_buses + 1):end] = ones(n_isolated_buses)
    A_J[(end - n_isolated_buses + 1):end] =
        [get_bus_index(x, bus_lookup, nr) for x in get_isolated_buses(ybus)]
    A_V[(end - n_isolated_buses + 1):end] = zeros(n_isolated_buses)
    data = SparseArrays.sparse(A_I, A_J, A_V)
    SparseArrays.dropzeros!(data)
    axes = (arc_ax, bus_ax)
    lookup = (make_ax_ref(arc_ax), make_ax_ref(bus_ax))
    subnetwork_axes = make_arc_bus_subnetwork_axes(ybus)

    return IncidenceMatrix(
        data,
        axes,
        lookup,
        subnetwork_axes,
        ybus.network_reduction_data,
    )
end

"""
Make subnetwork axes for LODF and VirtualLODF
"""
function make_arc_arc_subnetwork_axes(A::IncidenceMatrix)
    subnetwork_axes = Dict{Int, Tuple{Vector{Tuple{Int, Int}}, Vector{Tuple{Int, Int}}}}()
    for key in keys(A.subnetwork_axes)
        subnetwork_axes[key] = (A.subnetwork_axes[key][1], A.subnetwork_axes[key][1])
    end
    return subnetwork_axes
end
