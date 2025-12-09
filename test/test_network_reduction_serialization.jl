@testset failfast = true "Test serialization of PTDF with NetworkReductionData" begin
    # Test 1: Serialization with RadialReduction
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")

    # Create PTDF with radial reduction
    radial_reduction = RadialReduction()
    ptdf_radial = PTDF(sys5; network_reductions = NetworkReduction[radial_reduction])

    # Serialize and deserialize
    path = mktempdir()
    filename = joinpath(path, "ptdf_radial.h5")
    to_hdf5(ptdf_radial, filename)
    @test isfile(filename)

    ptdf_radial_loaded = PTDF(filename)

    # Check that network reduction data was preserved
    nrd_orig = get_network_reduction_data(ptdf_radial)
    nrd_loaded = get_network_reduction_data(ptdf_radial_loaded)

    # Check bus mappings
    @test nrd_loaded.bus_reduction_map == nrd_orig.bus_reduction_map
    @test nrd_loaded.reverse_bus_search_map == nrd_orig.reverse_bus_search_map
    @test nrd_loaded.irreducible_buses == nrd_orig.irreducible_buses
    @test nrd_loaded.removed_buses == nrd_orig.removed_buses
    @test nrd_loaded.removed_arcs == nrd_orig.removed_arcs

    # Check added maps
    @test nrd_loaded.added_branch_map == nrd_orig.added_branch_map
    @test nrd_loaded.added_admittance_map == nrd_orig.added_admittance_map

    # Check reduction container
    @test nrd_loaded.reductions == nrd_orig.reductions
    @test has_radial_reduction(nrd_loaded)
    @test !has_degree_two_reduction(nrd_loaded)
    @test !has_ward_reduction(nrd_loaded)

    # Check direct_branch_name_map
    @test nrd_loaded.direct_branch_name_map == nrd_orig.direct_branch_name_map

    # Test 2: Serialization with DegreeTwoReduction
    ptdf_d2 = PTDF(sys5; network_reductions = NetworkReduction[DegreeTwoReduction()])

    filename_d2 = joinpath(path, "ptdf_d2.h5")
    to_hdf5(ptdf_d2, filename_d2)
    ptdf_d2_loaded = PTDF(filename_d2)

    nrd_d2_orig = get_network_reduction_data(ptdf_d2)
    nrd_d2_loaded = get_network_reduction_data(ptdf_d2_loaded)

    @test nrd_d2_loaded.bus_reduction_map == nrd_d2_orig.bus_reduction_map
    @test nrd_d2_loaded.reverse_bus_search_map == nrd_d2_orig.reverse_bus_search_map
    @test nrd_d2_loaded.removed_buses == nrd_d2_orig.removed_buses
    @test nrd_d2_loaded.removed_arcs == nrd_d2_orig.removed_arcs
    @test has_degree_two_reduction(nrd_d2_loaded)

    # Test 3: Serialization with multiple reductions
    ptdf_multi = PTDF(
        sys5;
        network_reductions = NetworkReduction[RadialReduction(), DegreeTwoReduction()],
    )

    filename_multi = joinpath(path, "ptdf_multi.h5")
    to_hdf5(ptdf_multi, filename_multi)
    ptdf_multi_loaded = PTDF(filename_multi)

    nrd_multi_orig = get_network_reduction_data(ptdf_multi)
    nrd_multi_loaded = get_network_reduction_data(ptdf_multi_loaded)

    @test nrd_multi_loaded.bus_reduction_map == nrd_multi_orig.bus_reduction_map
    @test nrd_multi_loaded.reverse_bus_search_map == nrd_multi_orig.reverse_bus_search_map
    @test nrd_multi_loaded.removed_buses == nrd_multi_orig.removed_buses
    @test nrd_multi_loaded.removed_arcs == nrd_multi_orig.removed_arcs
    @test nrd_multi_loaded.reductions == nrd_multi_orig.reductions
    @test has_radial_reduction(nrd_multi_loaded)
    @test has_degree_two_reduction(nrd_multi_loaded)

    # Test 4: Serialization without reductions (backward compatibility)
    ptdf_no_red = PTDF(sys5)

    filename_no_red = joinpath(path, "ptdf_no_red.h5")
    to_hdf5(ptdf_no_red, filename_no_red)
    ptdf_no_red_loaded = PTDF(filename_no_red)

    @test isempty(get_network_reduction_data(ptdf_no_red_loaded))
    @test ptdf_no_red == ptdf_no_red_loaded

    # Test 5: Test compression works with network reduction data
    for compress in (true, false)
        filename_comp = joinpath(path, "ptdf_comp_$(compress).h5")
        to_hdf5(ptdf_radial, filename_comp; compress = compress)
        ptdf_comp_loaded = PTDF(filename_comp)

        @test get_network_reduction_data(ptdf_comp_loaded).bus_reduction_map ==
              get_network_reduction_data(ptdf_radial).bus_reduction_map
        @test get_network_reduction_data(ptdf_comp_loaded).reductions ==
              get_network_reduction_data(ptdf_radial).reductions
    end
end

@testset "Test serialization preserves network topology information" begin
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")

    # Create PTDF with reductions
    ptdf_orig = PTDF(sys5; network_reductions = NetworkReduction[RadialReduction()])

    # Serialize and deserialize
    path = mktempdir()
    filename = joinpath(path, "ptdf_topology.h5")
    to_hdf5(ptdf_orig, filename)
    ptdf_loaded = PTDF(filename)

    # Verify that topology information is preserved
    nrd_orig = get_network_reduction_data(ptdf_orig)
    nrd_loaded = get_network_reduction_data(ptdf_loaded)

    # Check that we can still query which buses were removed
    @test nrd_loaded.removed_buses == nrd_orig.removed_buses

    # Check that we can still query bus mappings
    for (bus, reduced_buses) in nrd_orig.bus_reduction_map
        @test haskey(nrd_loaded.bus_reduction_map, bus)
        @test nrd_loaded.bus_reduction_map[bus] == reduced_buses
    end

    # Check that reverse bus search map works
    for (reduced_bus, parent_bus) in nrd_orig.reverse_bus_search_map
        @test haskey(nrd_loaded.reverse_bus_search_map, reduced_bus)
        @test nrd_loaded.reverse_bus_search_map[reduced_bus] == parent_bus
    end
end
