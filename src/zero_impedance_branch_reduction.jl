"""
    ZeroImpedanceBranchReduction <: NetworkReduction

Internal reduction that merges electrically equivalent buses connected by zero-impedance branches.

Applied automatically as the first mandatory reduction step during [`Ybus`](@ref) construction.
Zero-impedance arcs are identified by scanning the assembled Ybus off-diagonal entries: an arc
(i, j) is zero-impedance when `real(Y[i,j]) == 0.0` and `imag(Y[i,j]) >= ZERO_IMPEDANCE_BRANCH_YBUS_SUSCEPTANCE_THRESHOLD`.
Transformer arcs are excluded — only non-transformer branches (Lines, MonitoredLines,
DiscreteControlledACBranch, etc.) are eligible for reduction, matching PSSE behavior.

The from-bus of the arc (as stored in the branch maps) is the surviving bus after the merge,
unless the to-bus is listed in `irreducible_buses`, in which case the surviving side is flipped.
Bus merging is applied by row/column summation of the Ybus data matrix, which automatically cancels
the series admittance of the zero-impedance arc while preserving shunt contributions from both buses.

# Fields
- `irreducible_buses::Set{Int}`: Bus numbers that must not be eliminated. When one side of a
  zero-impedance arc is irreducible the merge is still performed but the irreducible bus is
  kept as the surviving bus. When both sides are irreducible the merge is skipped with a warning.
"""
@kwdef struct ZeroImpedanceBranchReduction <: NetworkReduction
    irreducible_buses::Set{Int} = Set{Int}()
end

get_irreducible_buses(zir::ZeroImpedanceBranchReduction) = zir.irreducible_buses

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

_is_transformer(::PSY.TwoWindingTransformer) = true
_is_transformer(::PSY.ACTransmission) = false
_any_transformer(parallel_br::AbstractBranchesParallel) =
    any(_is_transformer(br) for br in parallel_br)

function _build_transformer_arc_set(nrd::NetworkReductionData)
    transformer_arcs = Set{Tuple{Int, Int}}()
    union!(transformer_arcs, keys(nrd.transformer3W_map))
    for (arc, branch) in nrd.direct_branch_map
        _is_transformer(branch) && push!(transformer_arcs, arc)
    end
    for (arc, parallel_br) in nrd.parallel_branch_map
        _any_transformer(parallel_br) && push!(transformer_arcs, arc)
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
        iszero(real(v)) || continue
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
        irreducible = get_irreducible_buses(reduction)
        from_irred = from_no ∈ irreducible
        to_irred = to_no ∈ irreducible
        if from_irred && to_irred
            @warn "Zero-impedance branch between two irreducible buses $from_no and $to_no; skipping merge."
            continue
        elseif to_irred
            # Flip so the irreducible bus survives.
            from_no, to_no = to_no, from_no
        end

        _update_bus_maps!(nr.reverse_bus_search_map, nr.bus_reduction_map, to_no, from_no)
        push!(nr.removed_arcs, arc_key)
    end
    nr.merged_bus_pairs = copy(nr.reverse_bus_search_map)
    return nr
end
