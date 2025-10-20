"""
LCC admittance matrix

# Arguments
- `data::SparseArrays.SparseMatrixCSC{ComplexF32, Int}`:
        The arc admittance matrix in the from-to direction
- `axes<:NTuple{2, Dict}`:
        Tuple containing two vectors (the first one showing the arc tuples,
        the second showing the buses numbers).
- `lookup<:NTuple{2, Dict}`:
        Tuple containing two dictionaries, the first mapping the arc tuples
        and the second the buses with their enumerated indexes.
- `direction::Symbol`:
        Direction of admittance (:FromTo or :ToFrom)
"""
struct LCCAdmittanceMatrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{ComplexF32}
    data::SparseArrays.SparseMatrixCSC{ComplexF32, Int}
    axes::Ax
    lookup::L
    direction::Symbol
end

get_axes(M::LCCAdmittanceMatrix) = M.axes
get_lookup(M::LCCAdmittanceMatrix) = M.lookup
get_arc_axis(M::LCCAdmittanceMatrix) = M.axes[1]
get_arc_lookup(M::LCCAdmittanceMatrix) = M.lookup[1]
get_bus_axis(M::LCCAdmittanceMatrix) = M.axes[2]
get_bus_lookup(M::LCCAdmittanceMatrix) = M.lookup[2]
get_direction(M::LCCAdmittanceMatrix) = M.direction

function LCCAdmittanceMatrix(lccs::Vector{PSY.TwoTerminalLCCLine}, direction::Symbol)
    lcc_arcs = [PSY.get_arc(x) for x in lccs if PSY.get_available(x) == true]
    arc_ax = [(lcc.from.number, lcc.to.number) for lcc in lcc_arcs]
    arc_lookup = Dict([x => ix for (ix, x) in enumerate(arc_ax)])
    bus_from = [lcc.from.number for lcc in lcc_arcs]
    bus_to = [lcc.to.number for lcc in lcc_arcs]
    bus_ax = sort(unique(vcat(bus_from, bus_to)))
    bus_lookup = Dict([x => ix for (ix, x) in enumerate(bus_ax)])
    n_bus = length(bus_ax)
    n_arc = length(arc_ax)
    data = SparseArrays.sparse([n_arc], [n_bus], ComplexF32(0.0))
    axes = (arc_ax, bus_ax)
    lookup = (arc_lookup, bus_lookup)
    return LCCAdmittanceMatrix(data, axes, lookup, direction)
end
