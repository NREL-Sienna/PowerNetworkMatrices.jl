"""
    ThreeWindingTransformerWinding{T<:PSY.ThreeWindingTransformer} <: PSY.ACTransmission

Internal object representing a single winding of a three-winding transformer. Do not export.

This structure is used internally to decompose three-winding transformers into individual
winding components for network matrix construction and analysis purposes.

# Fields
- `transformer::T`: The parent three-winding transformer object
- `winding_number::Int`: The winding number (1, 2, or 3) that this object represents

# Note
This is an internal object and should not be constructed directly by users or added to a system.
"""
struct ThreeWindingTransformerWinding{T <: PSY.ThreeWindingTransformer} <:
       PSY.ACTransmission
    transformer::T
    winding_number::Int
end

get_transformer(tw::ThreeWindingTransformerWinding) = tw.transformer
get_winding_number(tw::ThreeWindingTransformerWinding) = tw.winding_number

function get_name(three_wt_winding::ThreeWindingTransformerWinding)
    transformer = get_transformer(three_wt_winding)
    winding = get_winding_number(three_wt_winding)
    return PSY.get_name(transformer) * "_winding_$winding"
end

function get_series_susceptance(segment::ThreeWindingTransformerWinding)
    tfw = get_transformer(segment)
    winding_int = get_winding_number(segment)
    return PSY.get_series_susceptances(tfw)[winding_int]
end

"""
    get_equivalent_rating(tw::ThreeWindingTransformerWinding)

Get the rating for a specific winding of a three-winding transformer.
Returns the winding-specific rating if non-zero, otherwise returns the parent transformer rating.
"""
function get_equivalent_rating(tw::ThreeWindingTransformerWinding)
    tfw = get_transformer(tw)
    winding_num = get_winding_number(tw)

    winding_rating = if winding_num == 1
        PSY.get_rating_primary(tfw)
    elseif winding_num == 2
        PSY.get_rating_secondary(tfw)
    elseif winding_num == 3
        PSY.get_rating_tertiary(tfw)
    else
        throw(ArgumentError("Invalid winding number: $winding_num"))
    end
    if winding_rating != 0.0
        return winding_rating
    elseif isnothing(PSY.get_rating(tfw))
        return 0.0
    else
        return PSY.get_rating(tfw)
    end
end

"""
    get_equivalent_available(tw::ThreeWindingTransformerWinding)

Get the availability status for a specific winding of a three-winding transformer.
Returns the availability status of the parent transformer.
"""
function get_equivalent_available(tw::ThreeWindingTransformerWinding)
    tfw = get_transformer(tw)
    winding_num = get_winding_number(tw)

    winding_status = if winding_num == 1
        PSY.get_available_primary(tfw)
    elseif winding_num == 2
        PSY.get_available_secondary(tfw)
    elseif winding_num == 3
        PSY.get_available_tertiary(tfw)
    else
        throw(ArgumentError("Invalid winding number: $winding_num"))
    end
    return winding_status
end

function get_arc_tuple(tr::ThreeWindingTransformerWinding)
    t3W = get_transformer(tr)
    arc_number = get_winding_number(tr)
    if arc_number == 1
        return (
            PSY.get_number(PSY.get_from(PSY.get_primary_star_arc(t3W))),
            PSY.get_number(PSY.get_to(PSY.get_primary_star_arc(t3W))),
        )
    elseif arc_number == 2
        return (
            PSY.get_number(PSY.get_from(PSY.get_secondary_star_arc(t3W))),
            PSY.get_number(PSY.get_to(PSY.get_secondary_star_arc(t3W))),
        )
    elseif arc_number == 3
        return (
            PSY.get_number(PSY.get_from(PSY.get_tertiary_star_arc(t3W))),
            PSY.get_number(PSY.get_to(PSY.get_tertiary_star_arc(t3W))),
        )
    else
        throw(error("Three-winding transformer arc number must be 1, 2, or 3"))
    end
end

function add_to_map(device::ThreeWindingTransformerWinding, filters::Dict)
    return add_to_map(get_transformer(device), filters)
end
