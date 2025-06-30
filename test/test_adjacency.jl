@testset "Test connected components" begin
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

@testset "Test find subnetworks" begin
    n = 11
    buses = 100 .+ collect(1:n)
    edge_inds = [(1, 2), (2, 3), (3, 1), # cycle
        (4, 5), (6, 7), (8, 4), (8, 6), # two short chains that merge.
        # 9 is isolated.
        (10, 11)]
    edges = [(buses[i], buses[j]) for (i, j) in edge_inds]
    A = SparseArrays.sparse(I(n))
    for (i, j) in edge_inds
        A[i, j] = 1
        A[j, i] = 1
    end
    test_subnetworks = PNM.find_subnetworks(A, buses)
    expected =
        [Set(100 .+ (1:3)), Set(100 .+ (4:8)), Set(100 .+ (9:9)), Set(100 .+ (10:11))]
    @test length(values(test_subnetworks)) == length(expected)
    for (k, v) in test_subnetworks
        @test k in v
    end
    for k in expected
        @test k in values(test_subnetworks)
    end
end
