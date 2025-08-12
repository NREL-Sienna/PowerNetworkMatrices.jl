@testset "Test populate_branch_maps_by_type!" begin
    # This tests the function populate_branch_maps_by_type! by rebuilding the original branch maps and testing the round trip.
    sys = PSB.build_system(PSSEParsingTestSystems, "psse_14_network_reduction_test_system")
    ybus = Ybus(sys; network_reductions = NetworkReduction[DegreeTwoReduction()])
    nrd = ybus.network_reduction_data
    PNM.populate_branch_maps_by_type!(nrd)
    all_branch_maps_by_type = nrd.all_branch_maps_by_type
    nrd_rebuild = NetworkReductionData()
    for (map_key, v1) in all_branch_maps_by_type
        for (type, v2) in v1
            for (entry, v3) in v2
                map = getproperty(nrd_rebuild, Symbol(map_key))
                map[entry] = v3
            end
        end
    end
    for entry in [
        :direct_branch_map,
        :reverse_direct_branch_map,
        :parallel_branch_map,
        :reverse_parallel_branch_map,
        :series_branch_map,
        :reverse_series_branch_map,
        :transformer3W_map,
        :reverse_transformer3W_map,
    ]
        original_map = getproperty(nrd, entry)
        rebuilt_map = getproperty(nrd_rebuild, entry)
        @test original_map == rebuilt_map
    end
end
