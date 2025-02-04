struct NetworkReduction
    bus_reduction_map::Dict{Int, Set{Int}}
    reverse_bus_search_map::Dict{Int, Int}
    removed_branches::Set{String}
    retained_branches::Set{String}
    virtual_admittances::Vector{Tuple{Int64, Int64, ComplexF64}}
end

get_bus_reduction_map(rb::NetworkReduction) = rb.bus_reduction_map
get_reverse_bus_search_map(rb::NetworkReduction) = rb.reverse_bus_search_map
get_removed_branches(rb::NetworkReduction) = rb.removed_branches
get_retained_branches(rb::NetworkReduction) = rb.retained_branches
get_virtual_admittances(rb::NetworkReduction) = rb.virtual_admittances

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
    virtual_admittances::Vector{Tuple{Int64, Int64, ComplexF64}} = Vector{
        Tuple{Int64, Int64, ComplexF64},
    }())
    return NetworkReduction(
        bus_reduction_map,
        reverse_bus_search_map,
        removed_branches,
        retained_branches,
        virtual_admittances,
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
