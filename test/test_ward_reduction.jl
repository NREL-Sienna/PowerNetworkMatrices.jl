
@testset "Ward reduction" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    Ybus5 = Ybus(sys)

    all_buses = [get_number(b) for b in get_components(Bus, sys)]
    internal_buses = [5]
    boundary_buses = [1, 4]
    study_buses = union(internal_buses, boundary_buses)

    reduced_Ybus = WardReduction(Ybus5, sys, study_buses)
    reduced_Ybus
    @test reduced_Ybus isa ReducedYbus

    #test that external buses have been eliminated 
    external_buses = setdiff(all_buses, study_buses)
    for external_bus in external_buses
        @test (external_bus ∉ reduced_Ybus.axes[1])
        @test (external_bus ∉ reduced_Ybus.axes[2])
    end

    #test that internal entries are unchanged 
    for i in internal_buses
        for j in internal_buses
            @test Ybus5[i, j] == reduced_Ybus[i, j]
        end
    end

    #test size of reduced ybus 
    n_full = length(all_buses)
    n_reduced = length(study_buses)
    @test size(reduced_Ybus) == (n_reduced, n_reduced)
    @test size(Ybus5) == (n_full, n_full)

    #test result for different order of system buses
    reduced_Ybus2 = WardReduction(Ybus5, sys, [1, 4, 5])
    for i in study_buses
        for j in study_buses
            @test reduced_Ybus[i, j] == reduced_Ybus2[i, j]
        end
    end
end
