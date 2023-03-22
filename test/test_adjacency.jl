@testset "test connected components" begin
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    sys10 = get_10bus_test_system()

    M = AdjacencyMatrix(sys5)
    subnetworks = find_subnetworks(M)
    @test length(subnetworks) == 1

    M = AdjacencyMatrix(sys10; check_connectivity = false)
    subnetworks = find_subnetworks(M)
    @test length(subnetworks) == 2
end
