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
