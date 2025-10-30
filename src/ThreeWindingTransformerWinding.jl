struct ThreeWindingTransformerWinding{T<:PSY.ThreeWindingTransformer} <: PSY.ACTransmission
    transformer::T
    winding_number::Int
end

get_transformer(tw::ThreeWindingTransformerWinding) = tw.transformer
get_winding_number(tw::ThreeWindingTransformerWinding) = tw.winding_number
