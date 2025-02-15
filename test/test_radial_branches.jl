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
    rb = get_radial_reduction(sys)
    @test get_bus_reduction_map(rb)[7] == Set([61, 8])
    @test get_removed_branches(rb) == Set(["tl", "Trans4"])
    @test get_reduction_type(rb) == NetworkReductionTypes.RADIAL
end

@testset "Radial Branches Large" begin
    sys =
        build_system(MatpowerTestSystems, "matpower_ACTIVSg10k_sys")
    rb = get_radial_reduction(sys)
    for (k, v) in get_bus_reduction_map(rb)
        @test k ∉ v
    end
end

@testset "Check reference bus in Radial Branches" begin
    for name in ["matpower_ACTIVSg2000_sys", "matpower_ACTIVSg10k_sys"]
        sys = build_system(MatpowerTestSystems, name)
        a_mat = IncidenceMatrix(sys)
        rb = get_radial_reduction(sys)
        leaf_buses = Int64[]
        for i in keys(rb.bus_reduction_map)
            append!(leaf_buses, collect(rb.bus_reduction_map[i]))
        end
        @test all(a_mat.ref_bus_positions .∉ leaf_buses)
    end
end

function build_hvdc_with_small_island()
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")
    bus15 = ACBus(;
        number = 15,
        name = "Bus 15",
        bustype = ACBusTypes.REF,
        angle = 0.0,
        magnitude = 1.0,
        voltage_limits = (min = 0.9, max = 1.05),
        base_voltage = 69.0,
    )
    bus16 = ACBus(;
        number = 16,
        name = "Bus 16",
        bustype = ACBusTypes.PV,
        angle = 0.0,
        magnitude = 1.0,
        voltage_limits = (min = 0.9, max = 1.05),
        base_voltage = 69.0,
    )
    bus17 = ACBus(;
        number = 17,
        name = "Bus 17",
        bustype = ACBusTypes.PQ,
        angle = 0.0,
        magnitude = 1.0,
        voltage_limits = (min = 0.9, max = 1.05),
        base_voltage = 69.0,
    )
    add_component!(sys, bus15)
    add_component!(sys, bus16)
    add_component!(sys, bus17)

    line17 = Line(;
        name = "Line17",
        available = true,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = Arc(; from = bus15, to = bus16),
        r = 0.00281, # Per-unit
        x = 0.0281, # Per-unit
        b = (from = 0.00356, to = 0.00356), # Per-unit
        rating = 2.0, # Line rating of 200 MVA / System base of 100 MVA
        angle_limits = (min = -0.7, max = 0.7),
    )
    add_component!(sys, line17)
    line18 = Line(;
        name = "Line18",
        available = true,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = Arc(; from = bus15, to = bus17),
        r = 0.00281, # Per-unit
        x = 0.0281, # Per-unit
        b = (from = 0.00356, to = 0.00356), # Per-unit
        rating = 2.0, # Line rating of 200 MVA / System base of 100 MVA
        angle_limits = (min = -0.7, max = 0.7),
    )
    add_component!(sys, line18)
    load16 = PowerLoad(;
        name = "Bus16",
        available = true,
        bus = bus16,
        active_power = 0.0, # Per-unitized by device base_power
        reactive_power = 0.0, # Per-unitized by device base_power
        base_power = 10.0, # MVA
        max_active_power = 1.0, # 10 MW per-unitized by device base_power
        max_reactive_power = 0.0,
    )
    add_component!(sys, load16)
    gen17 = ThermalStandard(;
        name = "Bus17",
        available = true,
        status = true,
        bus = bus17,
        active_power = 0.0, # Per-unitized by device base_power
        reactive_power = 0.0, # Per-unitized by device base_power
        rating = 1.0, # 30 MW per-unitized by device base_power
        active_power_limits = (min = 0.2, max = 1.0), # 6 MW to 30 MW per-unitized by device base_power
        reactive_power_limits = nothing, # Per-unitized by device base_power
        ramp_limits = (up = 0.2, down = 0.2), # 6 MW/min up or down, per-unitized by device base_power
        operation_cost = ThermalGenerationCost(nothing),
        base_power = 30.0, # MVA
        time_limits = (up = 8.0, down = 8.0), # Hours
        must_run = false,
        prime_mover_type = PrimeMovers.CC,
        fuel = ThermalFuels.NATURAL_GAS,
    )
    add_component!(sys, gen17)
    bus14 = get_component(ACBus, sys, "Bus 14")
    hvdc1 = TwoTerminalHVDCLine(;
        name = "Line18",
        available = true,
        active_power_flow = 0.0,
        arc = Arc(; from = bus14, to = bus15),
        active_power_limits_from = (min = -100.0, max = 100.0),
        active_power_limits_to = (min = -100.0, max = 100.0),
        reactive_power_limits_from = (min = -100.0, max = 100.0),
        reactive_power_limits_to = (min = -100.0, max = 100.0),
    )
    add_component!(sys, hvdc1)
    return sys
end

@testset "Small island exception for radial reduction" begin
    sys = build_hvdc_with_small_island()
    rr = get_radial_reduction(sys)
    @test isa(PTDF(sys; network_reduction = rr), PTDF)
end
