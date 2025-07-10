@testset "Large 2d Reduction Test" begin
    sys =
        @test_logs (:error, r"no active generators found at bus") match_mode = :any build_system(
            MatpowerTestSystems,
            "matpower_ACTIVSg10k_sys",
        )

    A = AdjacencyMatrix(sys; check_connectivity = false)
    reduction = get_reduction(A, sys, Val(NetworkReductionTypes.TWO_D))
end
