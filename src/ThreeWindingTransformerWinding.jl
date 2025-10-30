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
struct ThreeWindingTransformerWinding{T<:PSY.ThreeWindingTransformer} <: PSY.ACTransmission
    transformer::T
    winding_number::Int
end

get_transformer(tw::ThreeWindingTransformerWinding) = tw.transformer
get_winding_number(tw::ThreeWindingTransformerWinding) = tw.winding_number
