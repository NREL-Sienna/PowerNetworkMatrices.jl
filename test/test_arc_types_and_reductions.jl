function make_3WT()
    # Create the system
    sys = System(100.0)

    # Create three buses
    bus1 = ACBus(
        1,
        "HV",
        true,
        "REF",
        0.0,
        1.0,
        (min = 0.9, max = 1.1),
        110.0,
        nothing,
        nothing,
    )
    bus2 =
        ACBus(2, "MV", true, "PQ", 0.0, 1.0, (min = 0.9, max = 1.1), 33.0, nothing, nothing)
    bus3 =
        ACBus(3, "LV", true, "PQ", 0.0, 1.0, (min = 0.9, max = 1.1), 11.0, nothing, nothing)

    add_component!(sys, bus1)
    add_component!(sys, bus2)
    add_component!(sys, bus3)

    # Create the star (hidden) bus for the transformer
    starbus = ACBus(
        999,
        "starbus",
        true,
        "PQ",
        0.0,
        1.0,
        (min = 0.9, max = 1.1),
        11.0,
        nothing,
        nothing,
    )
    add_component!(sys, starbus)

    # Create arcs from each bus to the star bus
    arc1 = Arc(bus1, starbus)
    arc2 = Arc(bus2, starbus)
    arc3 = Arc(bus3, starbus)

    # Create the three-winding transformer
    tf3w = Transformer3W(;
        name = "T1",
        available = true,
        primary_star_arc = arc1,
        secondary_star_arc = arc2,
        tertiary_star_arc = arc3,
        star_bus = starbus,
        active_power_flow_primary = 0.0,
        reactive_power_flow_primary = 0.0,
        active_power_flow_secondary = 0.0,
        reactive_power_flow_secondary = 0.0,
        active_power_flow_tertiary = 0.0,
        reactive_power_flow_tertiary = 0.0,
        r_primary = 0.01,
        x_primary = 0.1,
        r_secondary = 0.01,
        x_secondary = 0.1,
        r_tertiary = 0.01,
        x_tertiary = 0.1,
        rating = 100.0,
        r_12 = 0.01,
        x_12 = 0.1,
        r_23 = 0.01,
        x_23 = 0.1,
        r_13 = 0.01,
        x_13 = 0.1,
        base_power_12 = 100.0,
        base_power_23 = 100.0,
        base_power_13 = 100.0,
    )

    add_component!(sys, tf3w)
    return sys
end

"""A simple system where a 3WT bus-to-star arc is involved in a radial reduction."""
function make_3WT_radial()
    sys = make_3WT()

    hv = get_component(PSY.ACBus, sys, "HV")
    _add_simple_thermal_standard!(sys, hv, 0.5, 0.0)

    lv = get_component(PSY.ACBus, sys, "LV")
    low_voltage = get_base_voltage(lv)
    tail = _add_simple_bus!(sys, 4, "TAIL", ACBusTypes.PQ, low_voltage)

    _add_simple_line!(sys, lv, tail)
    _add_simple_load!(sys, tail, 1.0, 0.2)
    return sys
end

"""A simple system where a 3WT bus-to-star arc is involved in a degree-2 reduction."""
function make_3WT_deg2()
    sys = make_3WT()

    hv = get_component(PSY.ACBus, sys, "HV")
    lv = get_component(PSY.ACBus, sys, "LV")
    low_voltage = get_base_voltage(lv)

    last = _add_simple_bus!(sys, 4, "LAST", ACBusTypes.PQ, low_voltage)
    mid = _add_simple_bus!(sys, 5, "MID", ACBusTypes.PQ, low_voltage)
    _add_simple_line!(sys, lv, mid)
    _add_simple_line!(sys, mid, last)

    _add_simple_load!(sys, last, 1.0, 0.2)
    _add_simple_thermal_standard!(sys, hv, 0.5, 0.0)
    return sys
end

"""A simple system where parallel lines are involved in a radial reduction."""
function make_parallel_radial()
    sys = System(100.0)
    b1 = _add_simple_bus!(sys, 1, "bus_1", ACBusTypes.REF, 100.0)
    b2 = _add_simple_bus!(sys, 2, "bus_2", ACBusTypes.PQ, 100.0)
    b3 = _add_simple_bus!(sys, 3, "bus_3", ACBusTypes.PQ, 100.0)

    # add a circuit of lines between the 3 buses so they can't be reduced.
    _add_simple_line!(sys, b1, b2)
    _add_simple_line!(sys, b2, b3)
    _add_simple_line!(sys, b1, b3)

    # add radial tail onto b3 with parallel lines
    b4 = _add_simple_bus!(sys, 4, "bus_4", ACBusTypes.PQ, 100.0)
    b5 = _add_simple_bus!(sys, 5, "bus_5", ACBusTypes.PQ, 100.0)
    _add_simple_line!(sys, b3, b4)
    _add_simple_line!(sys, b3, b4)
    _add_simple_line!(sys, b4, b5)

    _add_simple_thermal_standard!(sys, b1, 0.5, 0.0)
    _add_simple_load!(sys, b3, 1.0, 0.2)
    _add_simple_load!(sys, b5, 1.0, 0.2)
    return sys
end

"""A simple system where parallel lines are involved in a degree-2 reduction."""
function make_parallel_deg2()
    sys = System(100.0)
    b1 = _add_simple_bus!(sys, 1, "bus_1", ACBusTypes.REF, 100.0)
    b2 = _add_simple_bus!(sys, 2, "bus_2", ACBusTypes.PQ, 100.0)
    b3 = _add_simple_bus!(sys, 3, "bus_3", ACBusTypes.PQ, 100.0)
    b4 = _add_simple_bus!(sys, 4, "bus_4", ACBusTypes.PQ, 100.0)

    # b2 has degree 2 (b1–b2–b3 path)
    _add_simple_line!(sys, b1, b2)
    _add_simple_line!(sys, b2, b3)

    # parallel arc between b3 and b4
    _add_simple_line!(sys, b3, b4)
    _add_simple_line!(sys, b3, b4)

    _add_simple_load!(sys, b4, 1.0, 0.2)
    _add_simple_thermal_standard!(sys, b1, 0.5, 0.0)
    return sys
end

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
    for br in PSY.get_components(PSY.ACBranch, sys)
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

@testset "3WT + radial" begin
    # test that it runs without error
    sys = make_3WT_radial()
    network_reductions = NetworkReduction[RadialReduction()]
    Y = Ybus(sys; network_reductions = network_reductions)
    @test Y isa Any
    # test that the 3WT arc was actually reduced
    trf = first(get_components(PSY.Transformer3W, sys))
    trf_arcs = Tuple{Int, Int}[PNM.get_arc_tuple((trf, i)) for i in 1:3]
    nrd = PNM.get_network_reduction_data(Y)
    @test any(arc in PNM.get_removed_arcs(nrd) for arc in trf_arcs) ||
          any(reverse(arc) in PNM.get_removed_arcs(nrd) for arc in trf_arcs)
end

@testset "3WT + degree-2" begin
    # test that it runs without error
    sys = make_3WT_deg2()
    network_reductions = NetworkReduction[DegreeTwoReduction()]
    Y = Ybus(sys; network_reductions = network_reductions)
    @test Y isa Any
    # test that the 3WT arc was actually reduced
    nrd = PNM.get_network_reduction_data(Y)
    @test Tuple{PSY.Transformer3W, Int} in types_in_series_reduction(nrd)
end

@testset "Parallel lines + radial" begin
    # test that it runs without error
    sys = make_parallel_radial()
    network_reductions = NetworkReduction[RadialReduction()]
    Y = Ybus(sys; network_reductions = network_reductions)
    @test Y isa Any
    # test that the parallel lines were reduced
    nrd = PNM.get_network_reduction_data(Y)
    parallel_arc = find_parallel_arc(sys)
    @test parallel_arc in PNM.get_removed_arcs(nrd)
end

@testset "Parallel lines + degree-2" begin
    # test that it runs without error
    sys = make_parallel_deg2()
    network_reductions = NetworkReduction[DegreeTwoReduction()]
    Y = Ybus(sys; network_reductions = network_reductions)
    @test Y isa Any
    # test that the parallel lines were reduced
    nrd = PNM.get_network_reduction_data(Y)
    @test Set{PSY.ACTransmission} in types_in_series_reduction(nrd)
end
