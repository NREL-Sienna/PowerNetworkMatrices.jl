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
            0.0,
            (from = 0.0, to = 0.0),
            100.0,
            (0.0, 0.0),
        ),
    )
    rb = RadialBranches(IncidenceMatrix(sys))
    @test get_bus_reduction_map(rb)[7] == Set([61, 8])
    @test get_radial_branches(rb) == Set(["tl", "Trans4"])
end

@testset "Radial Branches Large" begin
    sys =
        build_system(MatpowerTestSystems, "matpower_ACTIVSg10k_sys"; add_forecasts = false)
    rb = RadialBranches(IncidenceMatrix(sys))
    for (k, v) in get_bus_reduction_map(rb)
        @test k ∉ v
    end
end
