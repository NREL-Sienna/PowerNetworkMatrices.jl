@testset "Invalid reduction combinations" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    @test_throws IS.DataFormatError("RadialReduction is applied twice to the same system") Ybus(
        sys;
        network_reductions = NetworkReduction[RadialReduction(), RadialReduction()],
    )
    @test_throws IS.DataFormatError(
        "When applying both RadialReduction and DegreeTwoReduction, RadialReduction must be applied first",
    ) Ybus(
        sys;
        network_reductions = NetworkReduction[DegreeTwoReduction(), RadialReduction()],
    )
    @test_throws IS.DataFormatError(
        "RadialReduction reduction is applied after Ward reduction. Ward reduction must be applied last.",
    ) Ybus(
        sys;
        network_reductions = NetworkReduction[WardReduction([1, 2, 4]), RadialReduction()],
    )
end

@testset "14 bus; default reductions" begin
    sys = PSB.build_system(PSSEParsingTestSystems, "psse_14_network_reduction_test_system")
    A = IncidenceMatrix(sys)
    ybus = Ybus(sys)
    nrd = get_network_reduction_data(ybus)
    @test nrd.irreducible_buses == Set{Int}()
    @test length(keys(nrd.bus_reduction_map)) == 15
    @test nrd.bus_reduction_map[112] == Set([113])
    @test nrd.bus_reduction_map[104] == Set([105])
    @test nrd.reverse_bus_search_map == Dict{Int, Int}(105 => 104, 113 => 112)
    @test length(keys(nrd.direct_branch_map)) == 14
    @test length(keys(nrd.parallel_branch_map)) == 1
    @test length(keys(nrd.series_branch_map)) == 0
    @test length(keys(nrd.transformer3W_map)) == 6
    @test nrd.removed_buses == Set{Int}()
    @test nrd.removed_arcs == Set([(112, 113), (104, 105)])
    @test Set(keys(nrd.added_admittance_map)) == Set{Int}()
    @test Set(keys(nrd.added_branch_map)) == Set{Tuple{Int, Int}}()
    @test Set(A.axes[1]) == union(
        Set(keys(nrd.direct_branch_map)),
        Set(keys(nrd.parallel_branch_map)),
        Set(keys(nrd.series_branch_map)),
        Set(keys(nrd.transformer3W_map)),
    )
    @test Set(A.axes[2]) == Set(keys(nrd.bus_reduction_map))
end

@testset "14 bus; + radial reduction" begin
    sys = PSB.build_system(PSSEParsingTestSystems, "psse_14_network_reduction_test_system")
    A = IncidenceMatrix(sys; network_reductions = NetworkReduction[RadialReduction()])
    ybus = Ybus(sys; network_reductions = NetworkReduction[RadialReduction()])
    nrd = get_network_reduction_data(ybus)
    @test nrd.irreducible_buses == Set{Int}()
    @test length(keys(nrd.bus_reduction_map)) == 13
    @test nrd.bus_reduction_map[112] == Set([113])
    @test nrd.bus_reduction_map[104] == Set([105])
    @test nrd.bus_reduction_map[1001] == Set([107, 108])
    @test nrd.reverse_bus_search_map ==
          Dict(113 => 112, 105 => 104, 108 => 1001, 107 => 1001)
    @test length(keys(nrd.direct_branch_map)) == 13
    @test length(keys(nrd.parallel_branch_map)) == 1
    @test length(keys(nrd.series_branch_map)) == 0
    @test length(keys(nrd.transformer3W_map)) == 5
    @test nrd.removed_buses == Set{Int}()
    @test nrd.removed_arcs == Set([(107, 108), (107, 1001), (112, 113), (104, 105)])
    @test Set(keys(nrd.added_admittance_map)) == Set{Int}()
    @test Set(keys(nrd.added_branch_map)) == Set{Tuple{Int, Int}}()
    @test Set(A.axes[1]) == union(
        Set(keys(nrd.direct_branch_map)),
        Set(keys(nrd.parallel_branch_map)),
        Set(keys(nrd.series_branch_map)),
        Set(keys(nrd.transformer3W_map)),
    )
    @test Set(A.axes[2]) == Set(keys(nrd.bus_reduction_map))
end

@testset "14 bus; + degree two reduction" begin
    sys = PSB.build_system(PSSEParsingTestSystems, "psse_14_network_reduction_test_system")
    A = IncidenceMatrix(sys; network_reductions = NetworkReduction[DegreeTwoReduction()])
    ybus = Ybus(sys; network_reductions = NetworkReduction[DegreeTwoReduction()])
    nrd = get_network_reduction_data(ybus)
    @test nrd.irreducible_buses == Set{Int}()
    @test length(keys(nrd.bus_reduction_map)) == 14
    @test nrd.bus_reduction_map[112] == Set([113])
    @test nrd.bus_reduction_map[104] == Set([105])
    @test nrd.reverse_bus_search_map == Dict(113 => 112, 105 => 104)
    @test length(keys(nrd.direct_branch_map)) == 13
    @test length(keys(nrd.parallel_branch_map)) == 1
    @test length(keys(nrd.series_branch_map)) == 1
    @test length(keys(nrd.transformer3W_map)) == 5
    @test nrd.removed_buses == Set{Int}(107)
    @test nrd.removed_arcs == Set([(107, 108), (107, 1001), (112, 113), (104, 105)])
    @test Set(keys(nrd.added_admittance_map)) == Set{Int}()
    @test Set(keys(nrd.added_branch_map)) == Set{Tuple{Int, Int}}()
    @test Set(A.axes[1]) == union(
        Set(keys(nrd.direct_branch_map)),
        Set(keys(nrd.parallel_branch_map)),
        Set(keys(nrd.series_branch_map)),
        Set(keys(nrd.transformer3W_map)),
    )
    @test Set(A.axes[2]) == Set(keys(nrd.bus_reduction_map))
    ybus_full = Ybus(sys)
    @test isapprox(ybus[108, 1001]^-1, (ybus_full[108, 107])^-1 + ybus_full[107, 1001]^-1)
end

