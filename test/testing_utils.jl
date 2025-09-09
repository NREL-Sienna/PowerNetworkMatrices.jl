"""These functions are used to create simple components for the tests to have more compact code"""

"""
    _check_name(sys::System, name::String, component_type::DataType)
    Check if the name is unique in the system. If not, append a number to the name.
"""
function _check_name(sys::System, name::String, component_type::DataType)
    # Check if the name is unique
    check = true
    i = 1
    while check
        if has_component(sys, component_type, name)
            i += 1
            name = name * "_$i"
        else
            check = false
        end
    end
    return name
end

"""    
    _add_simple_bus!(sys::System, number::Int, bus_type::ACBusTypes, base_voltage::Number, voltage_magnitude::Float64=1.0, voltage_angle::Float64=0.0)
    Simplified function to create and add a bus to the system with the given parameters.
"""
function _add_simple_bus!(
    sys::System,
    number::Int,
    name::String,
    bus_type::ACBusTypes,
    base_voltage::Number,
    voltage_magnitude::Float64 = 1.0,
    voltage_angle::Float64 = 0.0,
)
    bus = ACBus(;
        number = number,
        name = name,
        available = true,
        bustype = bus_type,
        angle = voltage_angle,
        magnitude = voltage_magnitude,
        voltage_limits = (0.0, 2.0),
        base_voltage = Float64(base_voltage),
    )
    add_component!(sys, bus)
    return bus
end

"""    
    _add_simple_load!(sys::System, bus::ACBus, active_power::Number, reactive_power::Number)
    Simplified function to create and add a load to the system with the given parameters.
"""
function _add_simple_load!(
    sys::System,
    bus::ACBus,
    active_power::Number,
    reactive_power::Number,
)
    load = PowerLoad(;
        name = _check_name(sys, "load_$(get_number(bus))", PowerLoad),
        available = true,
        bus = bus,
        active_power = Float64(active_power), # Per-unitized by device base_power
        reactive_power = Float64(reactive_power), # Per-unitized by device base_power
        base_power = 1.0, # MVA
        max_active_power = 100.0, # 10 MW per-unitized by device base_power
        max_reactive_power = 100.0,
    )

    add_component!(sys, load)
    return load
end

"""    
    _add_simple_thermal_standard!(sys::System, bus::ACBus, active_power::Number=0.0, reactive_power::Number=0.0)
    Simplified function to create and add a thermal standard generator to the system with the given parameters.
"""
function _add_simple_thermal_standard!(
    sys::System,
    bus::ACBus,
    active_power::Number,
    reactive_power::Number,
)
    gen = ThermalStandard(;
        name = _check_name(sys, "thermal_standard_$(get_number(bus))", ThermalStandard),
        available = true,
        status = true,
        bus = bus,
        active_power = Float64(active_power),
        reactive_power = Float64(reactive_power),
        rating = 1.0,
        active_power_limits = (min = 0, max = 1),
        reactive_power_limits = (min = -1, max = 1),
        ramp_limits = nothing,
        operation_cost = ThermalGenerationCost(nothing),
        base_power = 100.0,
        time_limits = nothing,
        prime_mover_type = PrimeMovers.OT,
        fuel = ThermalFuels.OTHER,
        services = Device[],
        dynamic_injector = nothing,
        ext = Dict{String, Any}(),
    )
    add_component!(sys, gen)
    return gen
end

"""    
    _add_simple_line!(sys::System, bus1::ACBus, bus2::ACBus, r::Float64=1e-3, x::Float64=1e-3, b::Float64=0.0)
    Simplified function to create and add a line to the system with the given parameters.
"""
function _add_simple_line!(
    sys::System,
    bus1::ACBus,
    bus2::ACBus,
    r::Float64 = 1e-3,
    x::Float64 = 1e-3,
    b::Float64 = 0.0,
)
    line = Line(;
        name = _check_name(sys, "line_$(get_number(bus1))_$(get_number(bus2))", Line),
        available = true,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = Arc(; from = bus1, to = bus2),
        r = r,
        x = x,
        b = (from = b / 2, to = b / 2),
        rating = 1.0,
        angle_limits = (min = -pi / 2, max = pi / 2),
    )
    add_component!(sys, line)
    return line
end
