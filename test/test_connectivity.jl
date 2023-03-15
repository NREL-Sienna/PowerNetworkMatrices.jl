@testset "Test connected networks" begin
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
    remove_components!(sys, Line)
    @test (@test_logs (:warn, "Principal connected component does not contain:") match_mode =
        :any validate_connectivity(sys)) == false
    @test validate_connectivity(sys; connectivity_method = dfs_connectivity) ==
          false
end
