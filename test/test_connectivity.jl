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

@testset "Test connected components" begin
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    M = Ybus(sys5)
    subnetworks = find_subnetworks(M)
    @test length(subnetworks) == 1

    sys10 = PSB.build_system(PSISystems, "2Area 5 Bus System")
    M = Ybus(sys10)
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

@testset "Test matrices for connectivity corner cases" begin
    sys = PSB.build_system(PSISystems, "HVDC_TWO_RTO_RTS_5min_sys")
    ybus = Ybus(sys)
    A = IncidenceMatrix(ybus)
    ptdf = PTDF(sys)
    lodf = LODF(sys)
    vptdf = VirtualPTDF(sys)
    vlodf = VirtualLODF(sys)

    @test length(ybus.subnetwork_axes) == 2
    @test keys(ybus.subnetwork_axes) == keys(ptdf.subnetwork_axes) ==
          keys(lodf.subnetwork_axes) == keys(vptdf.subnetwork_axes) ==
          keys(vlodf.subnetwork_axes)
    for k in keys(ptdf.subnetwork_axes)
        @test iszero([ptdf[x, k] for x in PNM.get_arc_axis(ptdf)])
    end
    ref_bus_numbers = [
        get_number(x) for
        x in PSY.get_components(x -> get_bustype(x) == ACBusTypes.REF, ACBus, sys)
    ]
    for ref_bus in ref_bus_numbers
        @test ref_bus ∈ keys(ybus.subnetwork_axes)
    end
end

@testset "Small island corner cases" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")
    ybus_1 = Ybus(sys)
    ptdf_1 = PTDF(sys)
    lodf_1 = LODF(sys)
    vptdf_1 = VirtualPTDF(sys)
    vlodf_1 = VirtualLODF(sys)

    sys = build_hvdc_with_single_bus_island()
    ybus_2 = Ybus(sys)
    ptdf_2 = PTDF(sys)
    lodf_2 = LODF(sys)
    vptdf_2 = VirtualPTDF(sys)
    vlodf_2 = VirtualLODF(sys)

    sys = build_hvdc_with_small_island()
    ybus_3 = Ybus(
        sys;
        network_reductions = NetworkReduction[RadialReduction(;
            irreducible_buses = collect(1:14),
        )],
    )
    ptdf_3 = PTDF(
        sys;
        network_reductions = NetworkReduction[RadialReduction(;
            irreducible_buses = collect(1:14),
        )],
    )
    lodf_3 = LODF(
        sys;
        network_reductions = NetworkReduction[RadialReduction(;
            irreducible_buses = collect(1:14),
        )],
    )
    vptdf_3 = VirtualPTDF(
        sys;
        network_reductions = NetworkReduction[RadialReduction(;
            irreducible_buses = collect(1:14),
        )],
    )
    vlodf_3 = VirtualLODF(
        sys;
        network_reductions = NetworkReduction[RadialReduction(;
            irreducible_buses = collect(1:14),
        )],
    )

    for i in ptdf_1.axes[1], j in ptdf_1.axes[2]
        @test ptdf_1[j, i] == ptdf_2[j, i] == ptdf_3[j, i]
    end
    for i in lodf_1.axes[1], j in lodf_1.axes[2]
        @test lodf_1[i, j] == lodf_2[i, j] == lodf_3[i, j]
    end
    for i in vptdf_1.axes[1], j in vptdf_1.axes[2]
        @test vptdf_1[i, j] == vptdf_2[i, j] == vptdf_3[i, j]
    end
    for i in vlodf_1.axes[1], j in vlodf_1.axes[2]
        @test vlodf_1[i, j] == vlodf_2[i, j] == vlodf_3[i, j]
    end
end

@testset "Subnetwork algorithms" begin
    sys = build_hvdc_with_small_island()
    ybus = @test_logs (:info, r"Finding subnetworks via iterative union find") match_mode =
        :any Ybus(sys)
    ybus =
        @test_logs (:info, r"Finding subnetworks via depth first search") match_mode = :any Ybus(
            sys;
            subnetwork_algorithm = depth_first_search,
        )

    sub_1 = PNM.find_subnetworks(
        ybus.data,
        ybus.axes[1];
        subnetwork_algorithm = depth_first_search,
    )
    sub_2 = PNM.find_subnetworks(
        ybus.data,
        ybus.axes[1];
        subnetwork_algorithm = iterative_union_find,
    )
    @test sub_1 == sub_2
end

@testset "Goderya sparse vs GraphBLAS implementations" begin
    import Random
    import SparseArrays

    # Fixed seed for reproducibility
    Random.seed!(42)

    # Create a 100-bus system with 3 connected components of sizes 40, 30, 30
    # Buses are randomly assigned to components (not contiguous ranges)
    n = 100
    component_sizes = [40, 30, 30]

    # Shuffle bus indices and partition into components
    shuffled_buses = Random.shuffle(1:n)
    components = Vector{Int}[]
    offset = 0
    for size in component_sizes
        push!(components, shuffled_buses[(offset + 1):(offset + size)])
        offset += size
    end

    # Build adjacency matrix with random edges within each component
    I_idx = Int[]
    J_idx = Int[]
    V = ComplexF64[]

    for nodes in components
        # Add edges to form a connected component (path + extra edges)
        # First, create a path to ensure connectivity
        for i in 1:(length(nodes) - 1)
            push!(I_idx, nodes[i])
            push!(J_idx, nodes[i + 1])
            push!(V, -1.0 + 0.0im)
            push!(I_idx, nodes[i + 1])
            push!(J_idx, nodes[i])
            push!(V, -1.0 + 0.0im)
        end
        # Add some random extra edges within the component
        for _ in 1:div(length(nodes), 2)
            i, j = rand(nodes, 2)
            if i != j
                push!(I_idx, i)
                push!(J_idx, j)
                push!(V, -0.5 + 0.0im)
                push!(I_idx, j)
                push!(J_idx, i)
                push!(V, -0.5 + 0.0im)
            end
        end
    end

    # Add diagonal elements (self-admittance)
    for i in 1:n
        push!(I_idx, i)
        push!(J_idx, i)
        push!(V, 5.0 + 0.0im)
    end

    ybus = SparseArrays.sparse(I_idx, J_idx, V, n, n)

    # Test both implementations
    result_sparse = PNM._goderya_sparse(ybus)
    result_graphblas = PNM._goderya_graphblas(ybus)

    # Both should return the same row indices (order may differ)
    @test sort(result_sparse) == sort(result_graphblas)

    # Verify we detected the islands (not all n^2 entries filled)
    @test length(unique(result_sparse)) < n^2
end
