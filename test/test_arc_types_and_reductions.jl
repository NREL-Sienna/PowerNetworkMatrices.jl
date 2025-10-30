function types_in_series_reduction(nrd::PNM.NetworkReductionData)
    types = Set{DataType}()
    for segments in values(PNM.get_series_branch_map(nrd))
        for comp in segments
            push!(types, typeof(comp))
        end
    end
    return types
end

function find_parallel_arc(sys::System)
    arcs_seen = Set{Tuple{Int, Int}}()
    for br in PSY.get_components(
        x -> !(typeof(x) <: PSY.ThreeWindingTransformer),
        PSY.ACBranch,
        sys,
    )
        arc = PNM.get_arc_tuple(br)
        if arc in arcs_seen
            return arc
        else
            push!(arcs_seen, arc)
        end
    end
    error("No parallel arcs found")
    return (-1, -1)
end

function test_all_subtypes(sys::System, network_reductions)
    for T in subtypes(PNM.PowerNetworkMatrix)
        # arc admittance matrix constructor has different args.
        (T == PNM.Ybus || T == PNM.ArcAdmittanceMatrix) && continue
        M = T(sys; network_reductions = deepcopy(network_reductions))
        # test that it runs without error
        @test M isa T
    end
    # return Ybus so we can inspect the network reduction data.
    Y = Ybus(sys; network_reductions = deepcopy(network_reductions))
    @test Y isa PNM.Ybus
    return Y
end

@testset "3WT + radial" begin
    sys = build_system(PSB.PSITestSystems, "case10_radial_series_reductions")
    network_reductions = NetworkReduction[RadialReduction()]
    Y = test_all_subtypes(sys, network_reductions)
    # test that the 3WT arc was actually reduced
    trf = first(get_components(PSY.Transformer3W, sys))
    trf_arcs = Tuple{Int, Int}[PNM.get_arc_tuple((trf, i)) for i in 1:3]
    nrd = PNM.get_network_reduction_data(Y)
    @test any(arc in PNM.get_removed_arcs(nrd) for arc in trf_arcs) ||
          any(reverse(arc) in PNM.get_removed_arcs(nrd) for arc in trf_arcs)
end

@testset "3WT + degree-2" begin
    sys = build_system(PSB.PSITestSystems, "case10_radial_series_reductions")
    network_reductions = NetworkReduction[DegreeTwoReduction()]
    Y = test_all_subtypes(sys, network_reductions)
    # test that the 3WT arc was actually reduced
    nrd = PNM.get_network_reduction_data(Y)
    @test Tuple{PSY.Transformer3W, Int} in types_in_series_reduction(nrd)
end

@testset "Parallel lines + radial" begin
    sys = build_system(PSB.PSITestSystems, "case10_radial_series_reductions")
    network_reductions = NetworkReduction[RadialReduction()]
    Y = test_all_subtypes(sys, network_reductions)
    # test that the parallel lines were reduced
    nrd = PNM.get_network_reduction_data(Y)
    parallel_arc = find_parallel_arc(sys)
    @test parallel_arc in PNM.get_removed_arcs(nrd)
end

@testset "Parallel lines + degree-2" begin
    sys = build_system(PSB.PSITestSystems, "case10_radial_series_reductions")
    network_reductions = NetworkReduction[DegreeTwoReduction()]
    Y = test_all_subtypes(sys, network_reductions)
    # test that the parallel lines were reduced
    nrd = PNM.get_network_reduction_data(Y)
    @test Set{PSY.Line} in types_in_series_reduction(nrd)
end
