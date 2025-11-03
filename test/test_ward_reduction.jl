function _basic_test_ward_reduction(sys, study_buses)
    ybus = Ybus(sys; network_reductions = NetworkReduction[WardReduction(study_buses)])
    wr = get_network_reduction_data(ybus)
    @test isa(wr, NetworkReductionData)
    @test get_reductions(wr) ==
          PNM.ReductionContainer(; ward_reduction = WardReduction(study_buses))
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
    @test Set(PNM.get_bus_axis(PTDF_)) == Set(study_buses)
    @test Set(PNM.get_arc_axis(PTDF_)) ==
          union(Set(added_branch_arcs), Set(direct_arcs), Set(parallel_arcs))

    LODF_ = LODF(sys; network_reductions = NetworkReduction[WardReduction(study_buses)])
    @test Set(LODF_.axes[1]) ==
          union(Set(added_branch_arcs), Set(direct_arcs), Set(parallel_arcs))
    @test Set(LODF_.axes[2]) ==
          union(Set(added_branch_arcs), Set(direct_arcs), Set(parallel_arcs))

    vPTDF_ =
        VirtualPTDF(sys; network_reductions = NetworkReduction[WardReduction(study_buses)])
    @test Set(PNM.get_bus_axis(vPTDF_)) == Set(study_buses)
    @test Set(PNM.get_arc_axis(vPTDF_)) ==
          union(Set(added_branch_arcs), Set(direct_arcs), Set(parallel_arcs))

    vLODF_ =
        VirtualLODF(sys; network_reductions = NetworkReduction[WardReduction(study_buses)])
    @test Set(vLODF_.axes[1]) ==
          union(Set(added_branch_arcs), Set(direct_arcs), Set(parallel_arcs))
    @test Set(vLODF_.axes[2]) ==
          union(Set(added_branch_arcs), Set(direct_arcs), Set(parallel_arcs))
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

@testset "Test similarity of PTDF with Ward" begin
    sys = PSB.build_system(PSB.PSIDSystems, "3 Bus Inverter Base")
    ptdf_3 = PTDF(sys)
    ptdf_2 = @test_logs (
        :warn,
        r"Equivalent branch computed during Ward reduction is in parallel with existing system branch",
    ) match_mode = :any PTDF(
        sys;
        network_reductions = NetworkReduction[WardReduction([101, 102])],
    )
    existing_line_susceptance = PSY.get_series_susceptance(
        ptdf_2.network_reduction_data.direct_branch_map[(101, 102)],
    )
    ward_line_susceptance =
        1 / imag(1 / (ptdf_2.network_reduction_data.added_branch_map[(101, 102)]))
    ward_multiplier =
        existing_line_susceptance / (existing_line_susceptance + ward_line_susceptance)
    @test abs(ptdf_3[(101, 102), 102] - ward_multiplier * ptdf_2[(101, 102), 102]) < 0.007
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
    @test_throws IS.DataFormatError get_network_reduction_data(
        Ybus(
            sys;
            network_reductions = NetworkReduction[WardReduction([1, 2, 3, 4, 5, 100])],
        ),
    )
    @test_throws IS.DataFormatError get_network_reduction_data(
        Ybus(sys; network_reductions = NetworkReduction[WardReduction([2, 3, 4])]),
    )
end
