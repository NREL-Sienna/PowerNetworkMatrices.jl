#NOTE: hardcoded for testing system 
function get_reduction(
    A::IncidenceMatrix,
    sys::PSY.System,
    ::Val{NetworkReductionTypes.DEGREE_TWO},
)
    @error "DEGREE TWO HARDCODED REDUCTION FOR 14 BUS SYSTEM"
    br1 = PSY.get_component(PSY.Line, sys, "BUS 101-BUS 115-i_1")
    br2 = PSY.get_component(PSY.Line, sys, "BUS 115-BUS 102-i_1")
    return NetworkReduction(
        Dict{Int, Set{Int}}(),
        Dict{Int, Int}(),
        Dict{Tuple{Int, Int}, PSY.Branch}(),
        Dict{PSY.Branch, Tuple{Int, Int}}(),
        Dict{Tuple{Int, Int}, Set{PSY.Branch}}(),
        Dict{PSY.Branch, Tuple{Int, Int}}(),
        Dict{Tuple{Int, Int}, Set{PSY.Branch}}((101, 102) => Set([br1, br2])),          # Map from NEW arc to the series branches that comprise it.
        Dict{PSY.Branch, Tuple{Int, Int}}(br1 => (101, 102), br2 => (101, 102)),
        Dict{Tuple{Int, Int}, Tuple{PSY.ThreeWindingTransformer, Int}}(),
        Dict{Tuple{PSY.ThreeWindingTransformer, Int}, Tuple{Int, Int}}(),
        Set{Int}([115]),
        Set{Tuple{Int, Int}}([(101, 115), (115, 102)]),                                 # Set of OLD arcs to be removed 
        Vector{NetworkReductionTypes}([NetworkReductionTypes.DEGREE_TWO]),
    )
end
