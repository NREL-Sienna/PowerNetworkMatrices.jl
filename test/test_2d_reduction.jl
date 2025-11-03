@testset "Large 2d Reduction Test" begin
    sys = build_system(
        MatpowerTestSystems,
        "matpower_ACTIVSg10k_sys",
    )
    ybus = Ybus(sys)
    reduction = PNM.get_reduction(ybus, sys, DegreeTwoReduction())
    @test !isempty(reduction.series_branch_map)
end
