"""
    NetworkReductionData

Mutable struct containing all data mappings and metadata for network reduction operations.
This structure tracks how buses and branches are mapped, combined, or eliminated during
network reduction algorithms.

# Fields
- `irreducible_buses::Set{Int}`: Buses that cannot be reduced
- `bus_reduction_map::Dict{Int, Set{Int}}`: Maps retained buses to sets of eliminated buses  
- `reverse_bus_search_map::Dict{Int, Int}`: Maps eliminated buses to their parent buses
- `direct_branch_map::Dict{Tuple{Int, Int}, PSY.ACTransmission}`: One-to-one branch mappings
- `reverse_direct_branch_map::Dict{PSY.ACTransmission, Tuple{Int, Int}}`: Reverse direct mappings
- `parallel_branch_map::Dict{Tuple{Int, Int}, Set{PSY.ACTransmission}}`: Parallel branch combinations
- `reverse_parallel_branch_map::Dict{PSY.ACTransmission, Tuple{Int, Int}}`: Reverse parallel mappings  
- `series_branch_map::Dict{Tuple{Int, Int}, Vector{Any}}`: Series branch combinations
- `reverse_series_branch_map::Dict{Any, Tuple{Int, Int}}`: Reverse series mappings
- `transformer3W_map::Dict{Tuple{Int, Int}, Tuple{PSY.ThreeWindingTransformer, Int}}`: Three-winding transformer mappings
- `reverse_transformer3W_map::Dict{Tuple{PSY.ThreeWindingTransformer, Int}, Tuple{Int, Int}}`: Reverse transformer mappings
- `removed_buses::Set{Int}`: Set of buses eliminated from the network
- `removed_arcs::Set{Tuple{Int, Int}}`: Set of arcs eliminated from the network  
- `added_admittance_map::Dict{Int, Complex{Float32}}`: Admittances added to buses during reduction
- `added_branch_map::Dict{Tuple{Int, Int}, Complex{Float32}}`: New branches created during reduction
- `all_branch_maps_by_type::Dict{String, Any}`: Branch mappings organized by component type
- `reductions::ReductionContainer`: Container tracking applied reduction algorithms
"""
@kwdef mutable struct NetworkReductionData
    irreducible_buses::Set{Int} = Set{Int}() # Buses that are not reduced in the network reduction
    bus_reduction_map::Dict{Int, Set{Int}} = Dict{Int, Set{Int}}() # Maps reduced bus to the set of buses it was reduced to
    reverse_bus_search_map::Dict{Int, Int} = Dict{Int, Int}()
    direct_branch_map::Dict{Tuple{Int, Int}, PSY.ACTransmission} =
        Dict{Tuple{Int, Int}, PSY.ACTransmission}()
    reverse_direct_branch_map::Dict{PSY.ACTransmission, Tuple{Int, Int}} =
        Dict{PSY.ACTransmission, Tuple{Int, Int}}()
    parallel_branch_map::Dict{Tuple{Int, Int}, Set{PSY.ACTransmission}} =
        Dict{Tuple{Int, Int}, Set{PSY.ACTransmission}}()
    reverse_parallel_branch_map::Dict{PSY.ACTransmission, Tuple{Int, Int}} =
        Dict{PSY.ACTransmission, Tuple{Int, Int}}()
    series_branch_map::Dict{Tuple{Int, Int}, Vector{Any}} =
        Dict{Tuple{Int, Int}, Vector{Any}}()
    reverse_series_branch_map::Dict{Any, Tuple{Int, Int}} =
        Dict{Any, Tuple{Int, Int}}()
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
    added_admittance_map::Dict{Int, Complex{Float32}} = Dict{Int, Complex{Float32}}()
    added_branch_map::Dict{Tuple{Int, Int}, Complex{Float32}} =
        Dict{Tuple{Int, Int}, Complex{Float32}}()
    all_branch_maps_by_type::Dict{String, Any} = Dict{String, Any}()
    reductions::ReductionContainer = ReductionContainer()
end

function populate_branch_maps_by_type!(nrd::NetworkReductionData)
    all_branch_maps_by_type = Dict(
        "direct_branch_map" =>
            Dict{Type{<:PSY.ACTransmission}, Dict{Tuple{Int, Int}, PSY.ACTransmission}}(),
        "reverse_direct_branch_map" =>
            Dict{Type{<:PSY.ACTransmission}, Dict{PSY.ACTransmission, Tuple{Int, Int}}}(),
        "parallel_branch_map" =>
            Dict{
                Type{<:PSY.ACTransmission},
                Dict{Tuple{Int, Int}, Set{PSY.ACTransmission}},
            }(),
        "reverse_parallel_branch_map" =>
            Dict{Type{<:PSY.ACTransmission}, Dict{PSY.ACTransmission, Tuple{Int, Int}}}(),
        "series_branch_map" =>
            Dict{Type{<:PSY.ACTransmission}, Dict{Tuple{Int, Int}, Vector{Any}}}(),
        "reverse_series_branch_map" =>
            Dict{Type{<:PSY.ACTransmission}, Dict{Any, Tuple{Int, Int}}}(),
        "transformer3W_map" => Dict{
            Type{<:PSY.ThreeWindingTransformer},
            Dict{
                Tuple{Int, Int},
                Tuple{PSY.ThreeWindingTransformer, Int},
            },
        }(),
        "reverse_transformer3W_map" => Dict{
            Type{<:PSY.ThreeWindingTransformer},
            Dict{
                Tuple{PSY.ThreeWindingTransformer, Int},
                Tuple{Int, Int},
            },
        }())

    for (k, v) in nrd.direct_branch_map
        map_by_type = get!(
            all_branch_maps_by_type["direct_branch_map"],
            _get_segment_type(v),
            Dict{Tuple{Int, Int}, PSY.ACTransmission}(),
        )
        map_by_type[k] = v
    end
    for (k, v) in nrd.reverse_direct_branch_map
        map_by_type = get!(
            all_branch_maps_by_type["reverse_direct_branch_map"],
            _get_segment_type(k),
            Dict{PSY.ACTransmission, Tuple{Int, Int}}(),
        )
        map_by_type[k] = v
    end
    for (k, v) in nrd.parallel_branch_map
        map_by_type = get!(
            all_branch_maps_by_type["parallel_branch_map"],
            _get_segment_type(v),
            Dict{Tuple{Int, Int}, Set{PSY.ACTransmission}}(),
        )
        map_by_type[k] = v
    end
    for (k, v) in nrd.reverse_parallel_branch_map
        map_by_type = get!(
            all_branch_maps_by_type["reverse_parallel_branch_map"],
            _get_segment_type(k),
            Dict{PSY.ACTransmission, Tuple{Int, Int}}(),
        )
        map_by_type[k] = v
    end
    for (k, v) in nrd.series_branch_map
        #Repeated entry for each type in series chain
        for x in v
            map_by_type = get!(
                all_branch_maps_by_type["series_branch_map"],
                _get_segment_type(x),
                Dict{Tuple{Int, Int}, Vector{Any}}(),
            )
            map_by_type[k] = v
        end
    end
    for (k, v) in nrd.reverse_series_branch_map
        map_by_type = get!(
            all_branch_maps_by_type["reverse_series_branch_map"],
            _get_segment_type(k),
            Dict{Tuple{Int, Int}, Vector{Any}}(),
        )
        map_by_type[k] = v
    end
    for (k, v) in nrd.transformer3W_map
        map_by_type = get!(
            all_branch_maps_by_type["transformer3W_map"],
            _get_segment_type(v),
            Dict{Tuple{Int, Int}, Vector{Any}}(),
        )
        map_by_type[k] = v
    end
    for (k, v) in nrd.reverse_transformer3W_map
        map_by_type = get!(
            all_branch_maps_by_type["reverse_transformer3W_map"],
            _get_segment_type(k),
            Dict{Tuple{Int, Int}, Vector{Any}}(),
        )
        map_by_type[k] = v
    end

    nrd.all_branch_maps_by_type = all_branch_maps_by_type
    return
end

_get_segment_type(::T) where {T <: PSY.ACBranch} = T
_get_segment_type(x::Set{PSY.ACTransmission}) = typeof(first(x))
_get_segment_type(x::Tuple{PSY.ThreeWindingTransformer, Int}) = typeof(first(x))

get_irreducible_buses(rb::NetworkReductionData) = rb.irreducible_buses
"""
    get_bus_reduction_map(rb::NetworkReductionData)

Get the bus reduction map from NetworkReductionData.

# Arguments
- `rb::NetworkReductionData`: The network reduction data

# Returns
- `Dict{Int, Set{Int}}`: Dictionary mapping retained buses to sets of removed buses
"""
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
"""
    get_reductions(rb::NetworkReductionData)

Get the reduction container from NetworkReductionData.

# Arguments
- `rb::NetworkReductionData`: The network reduction data

# Returns
- `ReductionContainer`: Container with the applied network reductions
"""
get_reductions(rb::NetworkReductionData) = rb.reductions

has_radial_reduction(rb::NetworkReductionData) = has_radial_reduction(rb.reductions)
has_degree_two_reduction(rb::NetworkReductionData) = has_degree_two_reduction(rb.reductions)
has_ward_reduction(rb::NetworkReductionData) = has_ward_reduction(rb.reductions)

function Base.isempty(rb::NetworkReductionData)
    for field in fieldnames(NetworkReductionData)
        if !isempty(getfield(rb, field))
            return false
        end
    end
    return true
end

"""
   get_retained_branches_names(network_reduction_data::NetworkReductionData)

Gets the branch names that are retained after network reduction. This method only returns the
branch names from non-three winding transformer branches that have a one-to-one correspondence with
arcs after the reduction. This does not include parallel branches or branches that have been reduced as
part of a series chain of degree two nodes.

# Arguments
- `network_reduction_data::NetworkReductionData`

# Returns
- `Vector{String}`: Vector of the retained branch names.
"""
function get_retained_branches_names(network_reduction_data::NetworkReductionData)
    return [
        PSY.get_name(branch) for
        branch in keys(network_reduction_data.reverse_direct_branch_map)
    ]
end

"""
   get_ac_transmission_types(network_reduction_data::NetworkReductionData)

Gets the concrete types of all AC transmission branches included in an instance of NetworkReductionData

# Arguments
- `network_reduction_data::NetworkReductionData`

# Returns
- `Set{DataType}`: Vector of the retained branch types.
"""
function get_ac_transmission_types(network_reduction_data::NetworkReductionData)
    direct_types = Set(typeof.(keys(network_reduction_data.reverse_direct_branch_map)))
    parallel_types =
        Set(typeof.(keys(network_reduction_data.reverse_parallel_branch_map)))
    series_types = Set(typeof.(keys(network_reduction_data.reverse_series_branch_map)))
    transformer_3w_devices =
        Set(
            first(tuple) for tuple in keys(network_reduction_data.reverse_transformer3W_map)
        )
    transformer_3W_types = typeof.(transformer_3w_devices)
    return union(direct_types, parallel_types, series_types, transformer_3W_types)
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

"""
Interface to obtain the parent bus number of a reduced bus when radial branches are eliminated

# Arguments
- `rb::NetworkReduction`: NetworkReduction object
- `bus_number::Int`: Bus number of the reduced bus
"""
function get_mapped_bus_number(rb::NetworkReductionData, bus_number::Int)
    return get(rb.reverse_bus_search_map, bus_number, bus_number)
end

"""
Interface to obtain the parent bus number of a reduced bus when radial branches are eliminated

# Arguments
- `rb::NetworkReduction`: NetworkReduction object
- `bus::ACBus`: Reduced bus
"""
function get_mapped_bus_number(rb::NetworkReductionData, bus::PSY.ACBus)
    return get_mapped_bus_number(rb, PSY.get_number(bus))
end

"""
Interface to obtain the arc axis based on the network reduction data
"""
function get_arc_axis(nr::NetworkReductionData)
    direct_arcs = collect(keys(nr.direct_branch_map))
    parallel_arcs = collect(keys(nr.parallel_branch_map))
    series_arcs = collect(keys(nr.series_branch_map))
    transformer_arcs = collect(keys(nr.transformer3W_map))
    additional_arcs = collect(keys(nr.added_branch_map))
    arc_ax = unique(
        vcat(direct_arcs, parallel_arcs, series_arcs, transformer_arcs, additional_arcs),
    )
    return arc_ax
end

function is_arc_in_series_map(nr::NetworkReductionData, arc::Tuple{Int64, Int64})
    return haskey(nr.series_branch_map, arc)
end

function get_mapped_series_branch(nr::NetworkReductionData, arc::Tuple{Int64, Int64})
    if is_arc_in_series_map(nr, arc)
        return nr.series_branch_map[arc]
    else
        error("Arc $arc not found in series branch map")
    end
    return
end

function Base.show(io::IO, ::MIME{Symbol("text/plain")}, nrd::NetworkReductionData)
    println("Network Reduction Summary:")
    println("\tNumber of remapped buses: $(length(nrd.reverse_bus_search_map))")
    println("\tNumber of direct branch mappings: $(length(nrd.direct_branch_map))")
    println(
        "\tNumber of parallel arcs (number of branches): $(length(nrd.parallel_branch_map)) ($(length(nrd.reverse_parallel_branch_map)))",
    )
    println(
        "\tNumber of series arcs (number of branches): $(length(nrd.series_branch_map)) ($(length(nrd.reverse_series_branch_map)))",
    )
    println("\tNumber of 3WT winding arcs:$(length(nrd.transformer3W_map))")
    println("\tNumber of removed buses: $(length(nrd.removed_buses))")
    println("\tNumber of removed arcs: $(length(nrd.removed_arcs))")
    println("\tNumber of added branches: $(length(nrd.added_branch_map))")
    println("\tNumber of added admittances: $(length(nrd.added_admittance_map))")
end
