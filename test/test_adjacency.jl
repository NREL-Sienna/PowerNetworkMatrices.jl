@testset "test connected components" begin
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    M = AdjacencyMatrix(sys5)
    subnetworks = find_subnetworks(M)
    @test length(subnetworks) == 1

    sys10 = PSB.build_system(PSISystems, "2Area 5 Bus System")
    M = AdjacencyMatrix(sys10; check_connectivity = false)
    subnetworks_m = find_subnetworks(M)
    @test length(subnetworks_m) == 2
    @test all([6, 1] .∈ keys(subnetworks_m))

    subnetworks_sys = find_subnetworks(sys10)
    @test all([4, 9] .∈ keys(subnetworks_sys))
end
