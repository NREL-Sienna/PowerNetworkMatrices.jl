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

@testset "Basic degree two NetworkReductionData" begin
    sys = PSB.build_system(PSSEParsingTestSystems, "psse_14_network_reduction_test_system")
    ybus = Ybus(sys; check_connectivity = false)
    network_reduction_data = PNM.get_reduction(ybus, sys, DegreeTwoReduction())
    @test !isempty(PNM.get_series_branch_map(network_reduction_data))
    #@test !isempty(PNM.get_reverse_series_branch_map(network_reduction_data))  # FAILS - need to implement building reverse map
end

#TODO - add testing for new capabilities (i.e. degree two reduction)
sys = PSB.build_system(PSSEParsingTestSystems, "psse_14_network_reduction_test_system")
# Full Matrix 
ybus = Ybus(sys)
A = IncidenceMatrix(ybus)
nr = ybus.network_reduction_data;

# Radial Reduction
ybus_radial = Ybus(sys; network_reductions = NetworkReduction[RadialReduction()])
A_radial = IncidenceMatrix(ybus_radial)
nr_radial = ybus_radial.network_reduction_data;

# Degree Two Reduction 
#ybus_degree_two = Ybus(sys; network_reductions = NetworkReduction[DegreeTwoReduction()])
#A_degree_two = IncidenceMatrix(ybus_degree_two)
#nr_degree_two = ybus_degree_two.network_reduction_data;

# Radial + Degree Two Reduction 
#ybus_radial_degree_two = Ybus(
#    sys;
#    network_reductions = [RadialReduction(), DegreeTwoReduction()],
#)
#A_radial_degree_two = IncidenceMatrix(ybus_radial_degree_two)
#nr_radial_degree_two = ybus_radial_degree_two.network_reduction_data;

#@test length(nr_radial_degree_two.series_branch_map) == 1 #One set of series lines 
#@test length(nr_radial_degree_two.reverse_bus_search_map) == 4 #Two breaker/switches and a two radial arcs 
#@test length(nr_radial_degree_two.removed_arcs) == 6    #Two degree two connected arcs and the radial arc and the breaker/switches
#@test length(nr_radial_degree_two.removed_buses) == 1   #The degree two bus 
#@test length(nr_radial_degree_two.parallel_branch_map) == 1   #We have a parallel branch in original system

#Test modification of the ybus
#@test ybus[101, 102] == 0.0
#@test ybus_degree_two[101, 102] != 0.0
#@test ybus_degree_two[101, 102] == 1 / (1 / ybus[101, 115] + 1 / ybus[115, 102])
