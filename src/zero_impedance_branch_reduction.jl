"""
    ZeroImpedanceBranchReduction <: NetworkReduction

Internal reduction that merges electrically equivalent buses connected by zero-impedance branches.

Applied automatically as the first mandatory reduction step during [`Ybus`](@ref) construction.
Zero-impedance arcs are identified by scanning the assembled Ybus off-diagonal entries: an arc
(i, j) is zero-impedance when `real(Y[i,j]) == 0.0` and `imag(Y[i,j]) >= ZERO_IMPEDANCE_BRANCH_YBUS_SUSCEPTANCE_THRESHOLD`.
Transformer arcs are excluded — only non-transformer branches (Lines, MonitoredLines,
DiscreteControlledACBranch, etc.) are eligible for reduction, matching PSSE behavior.

The from-bus of the arc (as stored in the branch maps) is the surviving bus after the merge.
Bus merging is applied by row/column summation of the Ybus data matrix, which automatically cancels
the series admittance of the zero-impedance arc while preserving shunt contributions from both buses.
"""
@kwdef struct ZeroImpedanceBranchReduction <: NetworkReduction end

function _find_zero_impedance_arc_key(bus_i::Int, bus_j::Int, nrd::NetworkReductionData)
    for arc_tuple in ((bus_i, bus_j), (bus_j, bus_i))
        if haskey(nrd.direct_branch_map, arc_tuple) ||
           haskey(nrd.parallel_branch_map, arc_tuple) ||
           haskey(nrd.series_branch_map, arc_tuple)
            return arc_tuple
        end
    end
    return nothing
end

function _build_transformer_arc_set(nrd::NetworkReductionData)
    transformer_arcs = Set{Tuple{Int, Int}}()
    union!(transformer_arcs, keys(nrd.transformer3W_map))
    for (arc, branch) in nrd.direct_branch_map
        branch isa PSY.TwoWindingTransformer && push!(transformer_arcs, arc)
    end
    for (arc, parallel_br) in nrd.parallel_branch_map
        any(br isa PSY.TwoWindingTransformer for br in parallel_br) &&
            push!(transformer_arcs, arc)
    end
    return transformer_arcs
end

function get_reduction(
    ybus::Ybus,
    sys::PSY.System,
    reduction::ZeroImpedanceBranchReduction,
)
    nr = NetworkReductionData()
    nrd = get_network_reduction_data(ybus)
    bus_ax = get_bus_axis(ybus)
    transformer_arcs = _build_transformer_arc_set(nrd)
    nz_rows, nz_cols, nz_vals = SparseArrays.findnz(get_data(ybus))
    for (row_ix, col_ix, v) in zip(nz_rows, nz_cols, nz_vals)
        row_ix >= col_ix && continue
        real(v) == 0.0 || continue
        imag(v) >= ZERO_IMPEDANCE_BRANCH_YBUS_SUSCEPTANCE_THRESHOLD || continue
        bus_i = bus_ax[row_ix]
        bus_j = bus_ax[col_ix]
        # Transformer arcs are excluded from zero-impedance bus merging
        (bus_i, bus_j) ∈ transformer_arcs && continue
        (bus_j, bus_i) ∈ transformer_arcs && continue
        arc_key = _find_zero_impedance_arc_key(bus_i, bus_j, nrd)
        if arc_key === nothing
            @warn "Could not find branch map entry for zero-impedance arc between buses $bus_i and $bus_j; skipping."
            continue
        end
        from_no, to_no = arc_key

        _update_bus_maps!(nr.reverse_bus_search_map, nr.bus_reduction_map, to_no, from_no)
        push!(nr.removed_arcs, arc_key)
    end
    nr.merged_bus_pairs = copy(nr.reverse_bus_search_map)
    return nr
end
