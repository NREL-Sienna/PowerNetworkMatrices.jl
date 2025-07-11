@testset "Large 2d Reduction Test" begin
    sys =
        @test_logs (:error, r"no active generators found at bus") match_mode = :any build_system(
            MatpowerTestSystems,
            "matpower_ACTIVSg10k_sys",
        )
    ybus = Ybus(sys; check_connectivity = false)
    reduction = PNM.get_reduction(ybus, sys, Val(NetworkReductionTypes.DEGREE_TWO))
end
