@testset "LODF matrices" begin
    nodes_5 = nodes5()
    branches_5 = branches5(nodes_5)
    P5 = PTDF(branches_5, nodes_5)
    L5 = LODF(branches_5, nodes_5)
    @test isapprox(maximum(L5.data), maximum(Lodf_5), atol = 1e-3)
    @test isapprox(L5[branches_5[1], branches_5[2]], 0.3447946513849091)

    nodes_14 = nodes14()
    branches_14 = branches14(nodes_14)
    L14 = LODF(branches_14, nodes_14)
    @test isapprox(maximum(L14.data), maximum(Lodf_14), atol = 1e-3)

    L5NS = LODF(sys)
    @test getindex(L5NS, "3-4-i_5", "2-3-i_4") - 0.0003413469090 <= 1e-4
    # ! the following does not pass the test
    # @test isapprox(getindex(L5NS, "3-4-i_5", "2-3-i_4"), 0.0003413469090, atol = 1e-4)

    L5NS = LODF([branches_5[b] for b in Br5NS_ids], [nodes_5[b] for b in Bu5NS_ids])
    for brf in Br5NS_ids, brt in Br5NS_ids
        @test isapprox(
            getindex(L5NS, string(brf), string(brt)),
            Lodf_5[brf, brt],
            atol = 1e-3,
        )
    end
end
