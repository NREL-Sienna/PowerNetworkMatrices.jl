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
    @test isa(yft, BranchAdmittanceMatrix)
    @test isa(ytf, BranchAdmittanceMatrix)
    @test size(A) == size(yft) == size(ytf)
end
