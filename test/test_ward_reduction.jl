@testset "9bus ward reduction" begin
    sys = PSB.build_system(PSB.PSIDTestSystems, "psid_test_ieee_9bus")
    study_buses = [1, 2, 5, 4, 7]
    wr = get_ward_reduction(sys, study_buses)   #[4 and 7 = boundary; 3, 8, 9, 6  = external]
    external_buses =
        setdiff([get_number(x) for x in get_components(ACBus, sys)], study_buses)
    @test !isempty(wr.virtual_admittances)
    for external_bus in external_buses
        @test external_bus ∉ keys(wr.bus_reduction_map)
        @test external_bus ∈ keys(wr.reverse_bus_search_map)
    end
end

@testset "RTS ward reduction" begin
    sys = PSB.build_system(PSB.PSISystems, "RTS_GMLC_DA_sys")
    bus_numbers = [get_number(x) for x in get_components(ACBus, sys)]
    study_buses = filter!(x -> digits(x)[end] == 1, bus_numbers)  #study buses are from area 1 
    wr = get_ward_reduction(sys, study_buses)
    external_buses =
        setdiff([get_number(x) for x in get_components(ACBus, sys)], study_buses)
    @test !isempty(wr.virtual_admittances)
    for external_bus in external_buses
        @test external_bus ∉ keys(wr.bus_reduction_map)
        @test external_bus ∈ keys(wr.reverse_bus_search_map)
    end
end

