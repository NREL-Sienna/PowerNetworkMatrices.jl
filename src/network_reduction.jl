struct NetworkReduction
    bus_reduction_map::Dict{Int, Set{Int}}
    reverse_bus_search_map::Dict{Int, Int}
    removed_branches::Set{String}
    retained_branches::Set{String}
    added_branches::Vector{PSY.ACBranch}
    added_admittances::Vector{PSY.FixedAdmittance}
    reduction_type::Union{Nothing, Vector{NetworkReductionTypes}}
end

get_bus_reduction_map(rb::NetworkReduction) = rb.bus_reduction_map
get_reverse_bus_search_map(rb::NetworkReduction) = rb.reverse_bus_search_map
get_removed_branches(rb::NetworkReduction) = rb.removed_branches
get_retained_branches(rb::NetworkReduction) = rb.retained_branches
get_added_branches(rb::NetworkReduction) = rb.added_branches
get_added_admittances(rb::NetworkReduction) = rb.added_admittances
get_reduction_type(rb::NetworkReduction) = rb.reduction_type

function Base.isempty(rb::NetworkReduction)
    if !isempty(rb.bus_reduction_map)
        return false
    end
    if !isempty(rb.reverse_bus_search_map)
        return false
    end
    return true
end

function NetworkReduction(;
    bus_reduction_map::Dict{Int, Set{Int}} = Dict{Int, Set{Int}}(),
    reverse_bus_search_map::Dict{Int, Int} = Dict{Int, Int}(),
    removed_branches::Set{String} = Set{String}(),
    retained_branches::Set{String} = Set{String}(),
    added_branches::Vector{PSY.ACBranch} = Vector{PSY.ACBranch}(),
    added_admittances::Vector{PSY.FixedAdmittance} = Vector{PSY.FixedAdmittance}(),
    reduction_type::Union{Nothing, Vector{NetworkReductionTypes}} = nothing,
)
    return NetworkReduction(
        bus_reduction_map,
        reverse_bus_search_map,
        removed_branches,
        retained_branches,
        added_branches,
        added_admittances,
        reduction_type,
    )
end

function validate_reduction_type(
    reduction_type::NetworkReductionTypes,
    prior_reduction_types::Union{Nothing, Vector{NetworkReductionTypes}},
)
    if prior_reduction_types === nothing
        return
    else
        if reduction_type ∈ prior_reduction_types
            throw(IS.DataFormatError("$reduction_type is applied twice to the same system"))
        end
        if reduction_type == NetworkReductionTypes.BREAKER_SWITCH
            if NetworkReductionTypes.WARD ∈ prior_reduction_types
                throw(IS.DataFormatError("Cannot apply $reduction_type after Ward Reduction"))
            end
        elseif reduction_type == NetworkReductionTypes.RADIAL
            if NetworkReductionTypes.WARD ∈ prior_reduction_types
                throw(IS.DataFormatError("Cannot apply $reduction_type after Ward Reduction"))
            end
        elseif reduction_type == NetworkReductionTypes.WARD
        else
            error("Define validation for $reduction_type reduction")
        end
    end
end

function compose_reductions(nr_1::NetworkReduction, nr_2::NetworkReduction, n_buses::Int)
    reverse_map_2 = nr_2.reverse_bus_search_map
    bus_reduction_map = copy(nr_1.bus_reduction_map)
    for (k, v) in reverse_map_2
        if !haskey(bus_reduction_map, v)
            @warn "Bus $v was reduced in the current reduction and a previously applied reduction. Precedence is given to the prior reduction"
            continue
        else
            bus_reduction_map[v] = union(bus_reduction_map[v], k)
        end
    end
    reverse_bus_search_map = _make_reverse_bus_search_map(bus_reduction_map, n_buses)
    if nr_1.reduction_type === nothing
        reductions = nr_2.reduction_type
    else
        reductions = vcat(nr_1.reduction_type, nr_2.reduction_type)
    end
    return NetworkReduction(
        bus_reduction_map,
        reverse_bus_search_map,
        union(nr_1.removed_branches, nr_2.removed_branches),
        setdiff(nr_1.retained_branches, nr_2.removed_branches),
        vcat(nr_1.added_branches, nr_2.added_branches),
        vcat(nr_1.added_admittances, nr_2.added_admittances),
        reductions,
    )
end

##############################################################################
########################### Auxiliary functions ##############################
##############################################################################

function isequal(
    rb1::NetworkReduction,
    rb2::NetworkReduction,
)
    for field in fieldnames(typeof(rb1))
        if getfield(rb1, field) != getfield(rb2, field)
            return false
        end
    end
    return true
end
