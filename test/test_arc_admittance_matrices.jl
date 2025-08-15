@testset "First tests for branch admittance matrices" begin
    sys = PSB.build_system(PSSEParsingTestSystems, "psse_14_network_reduction_test_system")
    ybus = Ybus(
        sys;
        network_reductions = NetworkReduction[DegreeTwoReduction()],
        make_branch_admittance_matrices = true,
    )
    A = IncidenceMatrix(ybus)
    yft = ybus.branch_admittance_from_to
    ytf = ybus.branch_admittance_to_from
    @test isa(ybus, Ybus)
    @test isa(yft, ArcAdmittanceMatrix)
    @test isa(ytf, ArcAdmittanceMatrix)
    @test size(A) == size(yft) == size(ytf)

    sys = build_system(PSISystems, "RTS_GMLC_DA_sys")
    y = Ybus(sys; make_branch_admittance_matrices = true)
    Yft = y.branch_admittance_to_from
    @test length(PNM.get_arc_lookup(Yft)) == size(Yft.data)[1]
    @test length(PNM.get_bus_lookup(Yft)) == size(Yft.data)[2]
end
