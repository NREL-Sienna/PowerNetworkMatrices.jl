struct ThreeWindingTransformerWinding <: PSY.ACTransmission
    tranformer::PSY.ThreeWindingTransformer
    winding_number::Int
end

IS.@forward((ThreeWindingTransformerWinding, :tranformer), PSY.ThreeWindingTransformer)
get_winding_number(tw::ThreeWindingTransformerWinding) = tw.winding_number
