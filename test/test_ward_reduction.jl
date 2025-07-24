function _basic_test_ward_reduction(sys, study_buses)
    ybus = Ybus(sys; network_reductions = NetworkReduction[WardReduction(study_buses)])
    wr = get_network_reduction_data(ybus)
    @test isa(wr, NetworkReductionData)
    @test get_reductions(wr) == [WardReduction(study_buses)]
    external_buses =
        setdiff([get_number(x) for x in get_components(ACBus, sys)], study_buses)
    @test !isempty(wr.added_admittance_map)
    @test !isempty(wr.added_branch_map)
    for external_bus in external_buses
        @test external_bus ∉ keys(wr.bus_reduction_map)
        @test external_bus ∈ keys(wr.reverse_bus_search_map)
    end
end

function _test_matrices_ward_reduction(sys, study_buses)
    ybus = Ybus(sys; network_reductions = NetworkReduction[WardReduction(study_buses)])
    wr = get_network_reduction_data(ybus)
    added_branch_arcs = [x for x in keys(wr.added_branch_map)]
    direct_arcs = [x for x in keys(wr.direct_branch_map)]
    parallel_arcs = [x for x in keys(wr.parallel_branch_map)]
    A = IncidenceMatrix(
        sys;
        network_reductions = NetworkReduction[WardReduction(study_buses)],
    )
    @test Set(A.axes[1]) ==
          union(Set(added_branch_arcs), Set(direct_arcs), Set(parallel_arcs))
    @test Set(A.axes[2]) == Set(study_buses)

    Adj = AdjacencyMatrix(
        sys;
        network_reductions = NetworkReduction[WardReduction(study_buses)],
    )
    @test Set(Adj.axes[1]) == Set(study_buses)
    @test Set(Adj.axes[2]) == Set(study_buses)

    BA = BA_Matrix(sys; network_reductions = NetworkReduction[WardReduction(study_buses)])
    @test Set(BA.axes[2]) ==
          union(Set(added_branch_arcs), Set(direct_arcs), Set(parallel_arcs))
    @test Set(BA.axes[1]) == Set(study_buses)

    Y = Ybus(sys; network_reductions = NetworkReduction[WardReduction(study_buses)])
    @test Set(Y.axes[1]) == Set(study_buses)
    @test Set(Y.axes[2]) == Set(study_buses)

    PTDF_ = PTDF(sys; network_reductions = NetworkReduction[WardReduction(study_buses)])
    @test Set(PTDF_.axes[1]) == Set(study_buses)
    @test Set(PTDF_.axes[2]) ==
          union(Set(added_branch_arcs), Set(direct_arcs), Set(parallel_arcs))

    LODF_ = LODF(sys; network_reductions = NetworkReduction[WardReduction(study_buses)])
    @test Set(LODF_.axes[1]) ==
          union(Set(added_branch_arcs), Set(direct_arcs), Set(parallel_arcs))
    @test Set(LODF_.axes[2]) ==
          union(Set(added_branch_arcs), Set(direct_arcs), Set(parallel_arcs))

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

    sys = PSB.build_system(PSSEParsingTestSystems, "psse_14_network_reduction_test_system")
    study_buses = [101, 114, 110, 111]
    _basic_test_ward_reduction(sys, study_buses)
    _test_matrices_ward_reduction(sys, study_buses)
end

#TODO - develop more robust testing for correctness of PTDF
@testset "Test similarity of PTDF with Ward" begin
    sys = PSB.build_system(PSB.PSIDSystems, "3 Bus Inverter Base")
    ptdf_3 = PTDF(sys)
    wr = get_network_reduction_data(
        Ybus(sys; network_reductions = NetworkReduction[WardReduction([101, 102])]),
    )
    ptdf_2 = PTDF(sys; network_reduction = wr)
    @test abs(ptdf_3["BUS 1-BUS 2-i_1", 102] - ptdf_2["BUS 1-BUS 2-i_1", 102]) < 0.0025
end

@testset "Ward corner cases" begin
    sys = build_hvdc_with_small_island()
    wr = get_network_reduction_data(
        Ybus(sys; network_reductions = NetworkReduction[WardReduction([1, 2, 3])]),
    )

    @test isa(wr, NetworkReductionData)
    @test length(wr.bus_reduction_map) == 3

    wr = get_network_reduction_data(
        Ybus(sys; network_reductions = NetworkReduction[WardReduction([1])]),
    )
    @test isa(wr, NetworkReductionData)
    @test length(wr.bus_reduction_map) == 1
    @test length(wr.added_branch_map) == 0
    @test length(wr.added_admittance_map) == 1

    wr =
        @test_logs (:error, r"no boundary buses found") match_mode = :any get_network_reduction_data(
            Ybus(sys; network_reductions = NetworkReduction[WardReduction([15, 16, 17])]),
        )
    @test isa(wr, NetworkReductionData)
    @test length(wr.added_branch_map) == 0
    @test length(wr.added_admittance_map) == 0

    @test_throws IS.DataFormatError get_network_reduction_data(
        Ybus(
            sys;
            network_reductions = NetworkReduction[WardReduction([1, 2, 3, 4, 5, 17])],
        ),
    )
    @test_throws ErrorException get_network_reduction_data(
        Ybus(
            sys;
            network_reductions = NetworkReduction[WardReduction([1, 2, 3, 4, 5, 100])],
        ),
    )
    @test_throws IS.DataFormatError get_network_reduction_data(
        Ybus(sys; network_reductions = NetworkReduction[WardReduction([2, 3, 4])]),
    )
end
