@testset "Test connected networks" begin
    sys = PSB.build_system(PSB.MatpowerTestSystems, "matpower_case5_sys")
    @test validate_connectivity(sys)
    @test(
        @test_logs (
            :info,
            "Validating connectivity with depth first search (network traversal)",
        ) match_mode = :any validate_connectivity(
            sys,
            connectivity_method = dfs_connectivity,
        )
    )
    @test length(collect(find_connected_components(sys))[1]) == 5
end

@testset "Test disconnected networks" begin
    sys = PSB.build_system(PSB.MatpowerTestSystems, "matpower_case5_sys")
    remove_components!(sys, Line)
    @test (@test_logs (:warn, "The system contains islands") match_mode =
        :any validate_connectivity(sys)) == false
    @test validate_connectivity(sys; connectivity_method = dfs_connectivity) ==
          false
end
