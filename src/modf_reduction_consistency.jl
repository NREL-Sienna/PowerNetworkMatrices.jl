# Builds the "exception list" of buses that a VirtualMODF's outaged/monitored
# branches need kept through reduction (a reduced-away endpoint can't be queried
# as an arc, silently corrupting post-contingency rows) and injects it into the
# reductions.

# Both arc endpoints of a branch. `_add_arc_buses_to_irreducible!` dispatches on
# branch type, so three-winding transformers (three star arcs, no single
# `get_arc`) work here too. Non-branch components hit the `::PSY.Component` no-op.
function _accumulate_protected_buses!(buses::Vector{Int}, branch::PSY.ACTransmission)
    _add_arc_buses_to_irreducible!(buses, branch)
    return
end

function _accumulate_protected_buses!(::Vector{Int}, ::PSY.Component)
    return
end

# Buses of the components an outage declares monitored (a `Set{Base.UUID}` on the
# outage). A stale UUID is surfaced via warning, not silently skipped.
function _accumulate_monitored_buses!(
    buses::Vector{Int},
    sys::PSY.System,
    outage::PSY.Outage,
)
    for uuid in PSY.get_monitored_components(outage)
        local component
        try
            component = IS.get_component(sys, uuid)
        catch e
            e isa ArgumentError || rethrow()
            @warn "Outage monitored component UUID $uuid not found in system; " *
                  "cannot protect it from reduction."
            continue
        end
        _accumulate_protected_buses!(buses, component)
    end
    return
end

# Bus -> zero-impedance merge survivor, read from an unreduced `Ybus` (main's exact
# map, already flattened). Don't re-derive the susceptance threshold: `1/1e-4`
# underflows it.
function _zero_impedance_survivor_map(sys::PSY.System)
    return get_reverse_bus_search_map(get_network_reduction_data(Ybus(sys)))
end

"""
    _collect_protected_buses(sys) -> Set{Int}

Return the set of bus numbers that must NOT be reduced so that, for every
`PSY.Outage` supplemental attribute, both its outaged component(s) and the
components it declares as monitored (`get_monitored_components`) remain
expressible as arcs in the reduced network. This is the VirtualMODF
"exception list".

Buses are mapped through the mandatory zero-impedance merge: marking both
endpoints of a zero-impedance branch irreducible makes `ZeroImpedanceBranchReduction`
skip the merge, leaving a singular ABA. Protecting the survivor instead keeps the
real arcs queryable.
"""
function _collect_protected_buses(sys::PSY.System)
    buses = Int[]
    for outage in PSY.get_supplemental_attributes(PSY.Outage, sys)
        for component in PSY.get_associated_components(sys, outage)
            _accumulate_protected_buses!(buses, component)
        end
        _accumulate_monitored_buses!(buses, sys, outage)
    end
    isempty(buses) && return Set{Int}()
    zi_map = _zero_impedance_survivor_map(sys)
    return Set{Int}(get(zi_map, b, b) for b in buses)
end

# Sorted, de-duplicated union of existing irreducible/study buses with `protected`.
_merge_irreducible(existing, protected::Set{Int}) =
    sort!(collect(union(Set(existing), protected)))

# A copy of `reduction` with its protected-bus set extended. Empty `protected`
# returns the input unchanged (no allocation on the no-outage path).
function _augment_reduction(reduction::RadialReduction, protected::Set{Int})
    isempty(protected) && return reduction
    return RadialReduction(;
        irreducible_buses = _merge_irreducible(reduction.irreducible_buses, protected),
    )
end

function _augment_reduction(reduction::DegreeTwoReduction, protected::Set{Int})
    isempty(protected) && return reduction
    return DegreeTwoReduction(;
        irreducible_buses = _merge_irreducible(reduction.irreducible_buses, protected),
        reduce_reactive_power_injectors = reduction.reduce_reactive_power_injectors,
    )
end

# Ward retains `study_buses`; a protected external bus must join the study area
# so its incident branch is not equivalenced away.
function _augment_reduction(reduction::WardReduction, protected::Set{Int})
    isempty(protected) && return reduction
    return WardReduction(_merge_irreducible(reduction.study_buses, protected))
end

"""
    _adjust_reductions_for_protection(reductions, protected) -> Vector{NetworkReduction}

Return a new reduction vector in which each reduction's protected-bus set has
been extended with `protected`.
"""
function _adjust_reductions_for_protection(
    reductions::Vector{NetworkReduction},
    protected::Set{Int},
)
    return NetworkReduction[_augment_reduction(r, protected) for r in reductions]
end
