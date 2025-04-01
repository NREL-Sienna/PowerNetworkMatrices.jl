function _basic_test_ward_reduction(sys, study_buses)
    wr = get_ward_reduction(sys, study_buses)
    @test isa(wr, NetworkReduction)
    @test get_reduction_type(wr) == [NetworkReductionTypes.WARD]
        external_buses =
        setdiff([get_number(x) for x in get_components(ACBus, sys)], study_buses)
    @test !isempty(wr.added_admittances)
    @test !isempty(wr.added_branches)
    for external_bus in external_buses
        @test external_bus ∉ keys(wr.bus_reduction_map)
        @test external_bus ∈ keys(wr.reverse_bus_search_map)
    end
end

function _test_matrices_ward_reduction(sys, study_buses)
    wr = get_ward_reduction(sys, study_buses)
    added_branch_names = [get_name(x) for x in wr.added_branches]

    A = IncidenceMatrix(sys; network_reduction = wr)
    @test Set(A.axes[1]) == union(Set(added_branch_names), Set(wr.retained_branches))
    @test Set(A.axes[2]) == Set(study_buses)

    Adj = AdjacencyMatrix(sys; network_reduction = wr)
    @test Set(Adj.axes[1]) == Set(study_buses)
    @test Set(Adj.axes[2]) == Set(study_buses)

    BA = BA_Matrix(sys; network_reduction = wr)
    @test Set(BA.axes[2]) == union(Set(added_branch_names), Set(wr.retained_branches))
    @test Set(BA.axes[1]) == Set(study_buses)

    Y = Ybus(sys; network_reduction = wr)
    @test Set(Y.axes[1]) == Set(study_buses)
    @test Set(Y.axes[2]) == Set(study_buses)

    PTDF_ = PTDF(sys; network_reduction = wr)
    @test Set(PTDF_.axes[1]) == Set(study_buses)
    @test Set(PTDF_.axes[2]) == union(Set(added_branch_names), Set(wr.retained_branches))

    LODF_ = LODF(sys; network_reduction = wr)
    @test Set(LODF_.axes[1]) == union(Set(added_branch_names), Set(wr.retained_branches))
    @test Set(LODF_.axes[2]) == union(Set(added_branch_names), Set(wr.retained_branches))

    #TODO - add virtual PTDF/ virtual LODF
end

@testset "Basic ward reduction" begin
    sys = PSB.build_system(PSB.PSIDSystems, "3 Bus Inverter Base")
    study_buses = [101, 102]
    _basic_test_ward_reduction(sys, study_buses)
    _test_matrices_ward_reduction(sys, study_buses)

    sys = PSB.build_system(PSB.PSIDTestSystems, "psid_test_ieee_9bus")
    study_buses = [1, 2, 5, 4, 7]
    _basic_test_ward_reduction(sys, study_buses)
    _test_matrices_ward_reduction(sys, study_buses)

    sys = PSB.build_system(PSB.PSISystems, "RTS_GMLC_DA_sys")
    bus_numbers = [get_number(x) for x in get_components(ACBus, sys)]
    study_buses = filter!(x -> digits(x)[end] == 1, bus_numbers)  #study buses are from area 1 
    _basic_test_ward_reduction(sys, study_buses)
    _test_matrices_ward_reduction(sys, study_buses)
end

#TODO - develop more robust testing for correctness of PTDF
@testset "Test similarity of PTDF with Ward" begin
    sys = PSB.build_system(PSB.PSIDSystems, "3 Bus Inverter Base")
    ptdf_3 = PTDF(sys)
    wr = get_ward_reduction(sys, [101, 102])
    ptdf_2 = PTDF(sys; network_reduction = wr)
    @test abs(ptdf_3["BUS 1-BUS 2-i_1", 102] - ptdf_2["BUS 1-BUS 2-i_1", 102]) < 0.0025
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
        name = "Line19",
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

@testset "Ward corner cases" begin
    sys = build_hvdc_with_small_island()
    wr = get_ward_reduction(sys, [1, 2, 3])
    @test isa(wr, NetworkReduction)
    @test length(wr.bus_reduction_map) == 3

    wr = get_ward_reduction(sys, [1])
    @test isa(wr, NetworkReduction)
    @test length(wr.bus_reduction_map) == 1
    @test length(wr.added_branches) == 0
    @test length(wr.added_admittances) == 1

    #TODO - verify this is desired behavior for eliminating entire island. 
    wr = get_ward_reduction(sys, [15, 16, 17])
    @test isa(wr, NetworkReduction)
    @test length(wr.added_branches) == 0
    @test length(wr.added_admittances) == 0

    #TODO - throw the correct error type 
    @test_throws IS.DataFormatError get_ward_reduction(sys, [1, 2, 3, 4, 5, 17])
    @test_throws IS.DataFormatError get_ward_reduction(sys, [1, 2, 3, 4, 5, 100])
    @test_throws IS.DataFormatError get_ward_reduction(sys, [2, 3, 4])
end
