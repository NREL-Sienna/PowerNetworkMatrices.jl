# Consistency between network reductions and the components a VirtualMODF must
# keep queryable. An outaged or monitored branch whose endpoint buses get
# reduced away can no longer be expressed as an arc, which silently corrupts
# post-contingency PTDF rows. These helpers build the "exception list" of buses
# that must survive reduction, and inject it into the requested reductions.

# Accumulate the buses a component contributes to the protected set.
# ACTransmission branches contribute both arc endpoints; other component types
# (generators, loads, shunts) do not constrain the DC sensitivity arcs and are
# ignored for protection purposes.
function _accumulate_protected_buses!(buses::Vector{Int}, branch::PSY.ACTransmission)
    _add_arc_buses_to_irreducible!(buses, branch)
    return
end

# Three-winding transformers are `<: ACTransmission` but have no single `get_arc`
# (they expose three star arcs), so they need a dedicated method — the
# ACTransmission method above would throw `MethodError(get_arc, ::Transformer3W)`.
# The star-arc helper protects all four buses.
function _accumulate_protected_buses!(
    buses::Vector{Int},
    transformer::PSY.ThreeWindingTransformer,
)
    _add_arc_buses_to_irreducible!(buses, transformer)
    return
end

function _accumulate_protected_buses!(::Vector{Int}, ::PSY.Component)
    return
end

# Accumulate the buses of every component an outage declares as monitored. The
# monitored set lives on the outage as a `Set{Base.UUID}`; resolve each through
# the system. A stale UUID (component removed after the outage was attached) is
# surfaced, not silently skipped.
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

"""
    _collect_protected_buses(sys) -> Set{Int}

Return the set of bus numbers that must NOT be reduced so that, for every
`PSY.Outage` supplemental attribute, both its outaged component(s) and the
components it declares as monitored (`get_monitored_components`) remain
expressible as arcs in the reduced network. This is the VirtualMODF
"exception list".
"""
function _collect_protected_buses(sys::PSY.System)
    buses = Int[]
    for outage in PSY.get_supplemental_attributes(PSY.Outage, sys)
        for component in PSY.get_associated_components(sys, outage)
            _accumulate_protected_buses!(buses, component)
        end
        _accumulate_monitored_buses!(buses, sys, outage)
    end
    return Set{Int}(buses)
end

# Return a copy of `reduction` whose protected/irreducible bus set is extended
# with `protected`. Empty `protected` returns the input unchanged (identity), so
# the common no-outage path allocates nothing.
function _augment_reduction(reduction::RadialReduction, protected::Set{Int})
    isempty(protected) && return reduction
    merged = sort!(collect(union(Set(reduction.irreducible_buses), protected)))
    return RadialReduction(; irreducible_buses = merged)
end

function _augment_reduction(reduction::DegreeTwoReduction, protected::Set{Int})
    isempty(protected) && return reduction
    merged = sort!(collect(union(Set(reduction.irreducible_buses), protected)))
    return DegreeTwoReduction(;
        irreducible_buses = merged,
        reduce_reactive_power_injectors = reduction.reduce_reactive_power_injectors,
    )
end

# Ward retains `study_buses`; a protected external bus must join the study area
# so its incident branch is not equivalenced away.
function _augment_reduction(reduction::WardReduction, protected::Set{Int})
    isempty(protected) && return reduction
    merged = sort!(collect(union(Set(reduction.study_buses), protected)))
    return WardReduction(merged)
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
