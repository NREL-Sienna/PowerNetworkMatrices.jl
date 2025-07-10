@kwdef mutable struct NetworkReduction
    irreducible_buses::Set{Int} = Set{Int}() # Buses that are not reduced in the network reduction
    bus_reduction_map::Dict{Int, Set{Int}} = Dict{Int, Set{Int}}() # Maps reduced bus to the set of buses it was reduced to
    reverse_bus_search_map::Dict{Int, Int} = Dict{Int, Int}()
    direct_branch_map::Dict{Tuple{Int, Int}, PSY.Branch} = Dict{Tuple{Int, Int}, PSY.Branch}()
    reverse_direct_branch_map::Dict{PSY.Branch, Tuple{Int, Int}} = Dict{PSY.Branch, Tuple{Int, Int}}()
    parallel_branch_map::Dict{Tuple{Int, Int}, Set{PSY.Branch}} = Dict{Tuple{Int, Int}, Set{PSY.Branch}}()
    reverse_parallel_branch_map::Dict{PSY.Branch, Tuple{Int, Int}} = Dict{PSY.Branch, Tuple{Int, Int}}()
    series_branch_map::Dict{Tuple{Int, Int}, Set{PSY.Branch}} = Dict{Tuple{Int, Int}, Set{PSY.Branch}}()
    reverse_series_branch_map::Dict{PSY.Branch, Tuple{Int, Int}} = Dict{PSY.Branch, Tuple{Int, Int}}()
    transformer3W_map::Dict{
        Tuple{Int, Int},
        Tuple{PSY.ThreeWindingTransformer, Int},
    } = Dict{Tuple{Int, Int}, Tuple{PSY.ThreeWindingTransformer, Int}}()
    reverse_transformer3W_map::Dict{
        Tuple{PSY.ThreeWindingTransformer, Int},
        Tuple{Int, Int},
    } = Dict{Tuple{PSY.ThreeWindingTransformer, Int}, Tuple{Int, Int}}()
    removed_buses::Set{Int} = Set{Int}()
    removed_arcs::Set{Tuple{Int, Int}} = Set{Tuple{Int, Int}}()
    reduction_type::Vector{NetworkReductionTypes} = Vector{NetworkReductionTypes}()
end

get_irreducible_buses(rb::NetworkReduction) = rb.irreducible_buses
get_bus_reduction_map(rb::NetworkReduction) = rb.bus_reduction_map
get_reverse_bus_search_map(rb::NetworkReduction) = rb.reverse_bus_search_map
get_direct_branch_map(rb::NetworkReduction) = rb.direct_branch_map
get_reverse_direct_branch_map(rb::NetworkReduction) = rb.reverse_direct_branch_map
get_parallel_branch_map(rb::NetworkReduction) = rb.parallel_branch_map
get_reverse_parallel_branch_map(rb::NetworkReduction) = rb.reverse_parallel_branch_map
get_series_branch_map(rb::NetworkReduction) = rb.series_branch_map
get_reverse_series_branch_map(rb::NetworkReduction) = rb.reverse_series_branch_map
get_transformer3W_map(rb::NetworkReduction) = rb.transformer3W_map
get_reverse_transformer3W_map(rb::NetworkReduction) = rb.reverse_transformer3W_map
get_reduction_type(rb::NetworkReduction) = rb.reduction_type

function Base.isempty(rb::NetworkReduction)
    for field in fieldnames(NetworkReduction)
        if !isempty(getfield(rb, field))
            return false
        end
    end
    return true
end

function validate_reduction_type(
    reduction_type::NetworkReductionTypes,
    prior_reduction_types::Vector{NetworkReductionTypes},
)
    if length(prior_reduction_types) == 0
        return
    else
        if reduction_type ∈ prior_reduction_types
            throw(IS.DataFormatError("$reduction_type is applied twice to the same system"))
        end
        if reduction_type == NetworkReductionTypes.BREAKER_SWITCH
            if NetworkReductionTypes.WARD ∈ prior_reduction_types
                throw(
                    IS.DataFormatError("Cannot apply $reduction_type after Ward Reduction"),
                )
            end
        elseif reduction_type == NetworkReductionTypes.RADIAL
            if NetworkReductionTypes.WARD ∈ prior_reduction_types
                throw(
                    IS.DataFormatError("Cannot apply $reduction_type after Ward Reduction"),
                )
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
    return NetworkReduction(
        bus_reduction_map,
        reverse_bus_search_map,
        union(nr_1.removed_branches, nr_2.removed_branches),
        setdiff(nr_1.retained_branches, nr_2.removed_branches),
        vcat(nr_1.added_branches, nr_2.added_branches),
        vcat(nr_1.added_admittances, nr_2.added_admittances),
        vcat(nr_1.reduction_type, nr_2.reduction_type),
    )
end

##############################################################################
########################### Auxiliary functions ##############################
##############################################################################

function isequal(
    rb1::NetworkReduction,
    rb2::NetworkReduction,
)
    for field in fieldnames(NetworkReduction)
        if getfield(rb1, field) != getfield(rb2, field)
            return false
        end
    end
    return true
end
