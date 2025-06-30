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

    #TODO - fails because no boundary buses are found 
    #wr = get_ward_reduction(sys, [15, 16, 17])
    #@test isa(wr, NetworkReduction)
    #@test length(wr.added_branches) == 0
    #@test length(wr.added_admittances) == 0

    #TODO - throw the correct error type 
    @test_throws IS.DataFormatError get_ward_reduction(sys, [1, 2, 3, 4, 5, 17])
    @test_throws IS.DataFormatError get_ward_reduction(sys, [1, 2, 3, 4, 5, 100])
    @test_throws IS.DataFormatError get_ward_reduction(sys, [2, 3, 4])
end
