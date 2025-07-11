abstract type NetworkReduction end

@kwdef mutable struct NetworkReductionData
    irreducible_buses::Set{Int} = Set{Int}() # Buses that are not reduced in the network reduction
    bus_reduction_map::Dict{Int, Set{Int}} = Dict{Int, Set{Int}}() # Maps reduced bus to the set of buses it was reduced to
    reverse_bus_search_map::Dict{Int, Int} = Dict{Int, Int}()
    direct_branch_map::Dict{Tuple{Int, Int}, PSY.Branch} =
        Dict{Tuple{Int, Int}, PSY.Branch}()
    reverse_direct_branch_map::Dict{PSY.Branch, Tuple{Int, Int}} =
        Dict{PSY.Branch, Tuple{Int, Int}}()
    parallel_branch_map::Dict{Tuple{Int, Int}, Set{PSY.Branch}} =
        Dict{Tuple{Int, Int}, Set{PSY.Branch}}()
    reverse_parallel_branch_map::Dict{PSY.Branch, Tuple{Int, Int}} =
        Dict{PSY.Branch, Tuple{Int, Int}}()
    series_branch_map::Dict{Tuple{Int, Int}, Set{PSY.Branch}} =
        Dict{Tuple{Int, Int}, Set{PSY.Branch}}()
    reverse_series_branch_map::Dict{PSY.Branch, Tuple{Int, Int}} =
        Dict{PSY.Branch, Tuple{Int, Int}}()
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
    added_admittance_map::Dict{Int, PSY.FixedAdmittance} = Dict{Int, PSY.FixedAdmittance}()
    added_branch_map::Dict{Tuple{Int, Int}, PSY.Line} = Dict{Tuple{Int, Int}, PSY.Line}()
    reductions::Vector{NetworkReduction} = Vector{NetworkReduction}() #store the reductions that were applied.
end

get_irreducible_buses(rb::NetworkReductionData) = rb.irreducible_buses
get_bus_reduction_map(rb::NetworkReductionData) = rb.bus_reduction_map
get_reverse_bus_search_map(rb::NetworkReductionData) = rb.reverse_bus_search_map
get_direct_branch_map(rb::NetworkReductionData) = rb.direct_branch_map
get_reverse_direct_branch_map(rb::NetworkReductionData) = rb.reverse_direct_branch_map
get_parallel_branch_map(rb::NetworkReductionData) = rb.parallel_branch_map
get_reverse_parallel_branch_map(rb::NetworkReductionData) = rb.reverse_parallel_branch_map
get_series_branch_map(rb::NetworkReductionData) = rb.series_branch_map
get_reverse_series_branch_map(rb::NetworkReductionData) = rb.reverse_series_branch_map
get_transformer3W_map(rb::NetworkReductionData) = rb.transformer3W_map
get_reverse_transformer3W_map(rb::NetworkReductionData) = rb.reverse_transformer3W_map
get_removed_buses(rb::NetworkReductionData) = rb.removed_buses
get_removed_arcs(rb::NetworkReductionData) = rb.removed_arcs
get_added_admittance_map(rb::NetworkReductionData) = rb.added_admittance_map
get_added_branch_map(rb::NetworkReductionData) = rb.added_branch_map
get_reductions(rb::NetworkReductionData) = rb.reductions

function Base.isempty(rb::NetworkReductionData)
    for field in fieldnames(NetworkReductionData)
        if !isempty(getfield(rb, field))
            return false
        end
    end
    return true
end

#TODO - refactor and use this function to validate that the order to reductions passed is a valid combination. 
#= function validate_reduction_type(
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
 =#

##############################################################################
########################### Auxiliary functions ##############################
##############################################################################

function isequal(
    rb1::NetworkReductionData,
    rb2::NetworkReductionData,
)
    for field in fieldnames(NetworkReductionData)
        if getfield(rb1, field) != getfield(rb2, field)
            return false
        end
    end
    return true
end

function Base.:(==)(x::T1, y::T1) where {T1 <: NetworkReduction}
    for field in fieldnames(T1)
        if getfield(x, field) != getfield(y, field)
            return false
        end
    end
    return true
end
