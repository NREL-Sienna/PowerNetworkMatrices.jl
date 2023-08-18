@testset "Test connected networks" begin
    sys = PSB.build_system(PSB.MatpowerTestSystems, "matpower_case5_sys")
    @test validate_connectivity(sys)
    @test(
        @test_logs (
            :info,
            "Validating connectivity with depth first search (network traversal)",
        ) match_mode = :any validate_connectivity(sys)
    )
end

@testset "Test disconnected networks" begin
    sys = PSB.build_system(PSB.MatpowerTestSystems, "matpower_case5_sys")
    remove_components!(sys, Line)
    @test(
        @test_logs (
            :warn,
            "Bus 1 is islanded",
        ) match_mode = :any validate_connectivity(sys) == false
    )
end