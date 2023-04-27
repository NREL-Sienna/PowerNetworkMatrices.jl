@testset "LODF matrices" begin
    """
    NOTE: LODF is transposed
    """

    # get system
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    nodes_5 = nodes5()
    branches_5 = branches5(nodes_5)

    # get LODF
    L5 = LODF(branches_5, nodes_5)
    @test isapprox(maximum(L5.data), maximum(Lodf_5), atol = 1e-3)
    @test isapprox(L5[branches_5[2], branches_5[1]], 0.3447946513849093)

    nodes_14 = nodes14()
    branches_14 = branches14(nodes_14)
    L14 = LODF(branches_14, nodes_14)
    @test isapprox(maximum(L14.data), maximum(Lodf_14), atol = 1e-3)

    L5NS = LODF(sys5)
    @test getindex(L5NS, "6", "5") - -0.3071 <= 1e-4
    total_error = abs.(L5NS.data' .- Lodf_5)
    @test isapprox(sum(total_error), 0.0, atol = 1e-3)

    A = IncidenceMatrix(sys5)
    P5 = PTDF(sys5)
    L5NS_from_ptdf = LODF(A, P5)
    @test getindex(L5NS_from_ptdf, "6", "5") - -0.3071 <= 1e-4
    total_error = abs.(L5NS_from_ptdf.data' .- Lodf_5)
    @test isapprox(sum(total_error), 0.0, atol = 1e-3)
end
