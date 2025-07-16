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
    series_branch_map::Dict{Tuple{Int, Int}, Vector{Any}} =
        Dict{Tuple{Int, Int}, Vector{Any}}()
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

function validate_reduction_type(
    reduction::T,
    prior_reductions::Vector{NetworkReduction},
) where {T <: NetworkReduction}
    prior_reduction_types = [typeof(x) for x in prior_reductions]
    if length(prior_reductions) == 0
        return
    else
        if T ∈ prior_reduction_types
            throw(IS.DataFormatError("$T is applied twice to the same system"))
        end
        if WardReduction ∈ prior_reduction_types
            throw(
                IS.DataFormatError(
                    "$T reduction is applied after Ward reduction. Ward reduction must be applied last.",
                ),
            )
        end
        if T == RadialReduction
            if DegreeTwoReduction ∈ prior_reduction_types
                throw(
                    IS.DataFormatError(
                        "When applying both RadialReduction and DegreeTwoReduction, RadialReduction must be applied first",
                    ),
                )
            end
        end
    end
end

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
