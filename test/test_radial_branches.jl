@testset "Radial Branches" begin
    sys = build_system(PSITestSystems, "c_sys14"; add_forecasts = false)
    n = first(get_components(ACBus, sys))
    n2 = deepcopy(n)
    n2.internal = PowerSystems.IS.InfrastructureSystemsInternal()
    set_name!(n2, "TestBus")
    set_number!(n2, 61)
    set_base_voltage!(n2, 18.0)
    add_component!(sys, n2)
    arc = Arc(get_component(ACBus, sys, "Bus 8"), n2)
    add_component!(sys, arc)
    add_component!(
        sys,
        Line(
            "tl",
            true,
            0.0,
            0.0,
            arc,
            0.0,
            0.01,    #cannot have zero impedance line
            (from = 0.0, to = 0.0),
            100.0,
            (0.0, 0.0),
        ),
    )
    Y = Ybus(sys; network_reductions = [NetworkReductionTypes.RADIAL])
    rb = get_network_reduction(Y)
    @test rb.bus_reduction_map[7] == Set([61, 8])
    @test rb.reverse_bus_search_map[61] == rb.reverse_bus_search_map[8] == 7
    @test length(rb.direct_branch_map) == 19
    @test length(rb.reverse_direct_branch_map) == 19
    @test length(rb.parallel_branch_map) == 0
    @test length(rb.reverse_parallel_branch_map) == 0
    @test length(rb.series_branch_map) == 0
    @test length(rb.reverse_series_branch_map) == 0
    @test length(rb.transformer3W_map) == 0
    @test length(rb.reverse_transformer3W_map) == 0
    @test length(rb.removed_buses) == 0
    @test rb.removed_arcs == Set([(7, 8), (8, 61)])
    @test get_reduction_type(rb) == [NetworkReductionTypes.RADIAL]
end

@testset "Radial Branches Large" begin
    sys =
        @test_logs (:error, r"no active generators found at bus") match_mode = :any build_system(
            MatpowerTestSystems,
            "matpower_ACTIVSg10k_sys",
        )
    Y = Ybus(sys; network_reductions = [NetworkReductionTypes.RADIAL])
    rb = get_network_reduction(Y)
    for (k, v) in get_bus_reduction_map(rb)
        @test k ∉ v
    end
end

@testset "Check reference bus in Radial Branches" begin
    for name in ["matpower_ACTIVSg2000_sys", "matpower_ACTIVSg10k_sys"]
        sys =
            @test_logs (:error, r"no active generators found at bus") match_mode = :any build_system(
                MatpowerTestSystems,
                name,
            )
        a_mat = IncidenceMatrix(sys)
        Y = Ybus(sys; network_reductions = [NetworkReductionTypes.RADIAL])
        rb = get_network_reduction(Y)
        leaf_buses = Int[]
        for i in keys(rb.bus_reduction_map)
            append!(leaf_buses, collect(rb.bus_reduction_map[i]))
        end
        leaf_positions = [a_mat.lookup[2][x] for x in leaf_buses]
        @test all(a_mat.ref_bus_positions .∉ leaf_positions)
    end
end

@testset "Small island exception for radial reduction" begin
    sys = build_hvdc_with_small_island()
    a_mat = IncidenceMatrix(sys; check_connectivity = false)
    rr = get_radial_reduction(a_mat)
    @test haskey(rr.reverse_bus_search_map, 16)
    @test haskey(rr.reverse_bus_search_map, 17)
    rr = get_radial_reduction(a_mat; exempt_buses = [16, 17])
    @test !haskey(rr.reverse_bus_search_map, 16)
    @test !haskey(rr.reverse_bus_search_map, 17)
    #@test isa(PTDF(sys; network_reduction = rr), PTDF)     #TODO - look at this test once PTDF works; I think I have captured the functionality with other added tests
end
