@testset failfast = true "Test VirtualLODF with radial lines" begin
    # get the system
    sys = PSB.build_system(PSB.PSITestSystems, "test_RTS_GMLC_sys")
    # get the PTDF matrix for reference
    lodf_rad = LODF(sys; reduce_radial_branches = true)
    vlodf_rad = VirtualLODF(sys; reduce_radial_branches = true)

    for i in axes(lodf_rad, 2)
        virtual = vlodf_rad[i, :]
        for j in axes(lodf_rad, 1)
            # check values using PTDFs axes
            @test isapprox(lodf_rad[i, j], vlodf_rad[i, j]; atol = 1e-10)
        end
    end
    # Check the cache is populated
    @test length(vlodf_rad.cache) == length(vlodf_rad.axes[1])
    @test length(vlodf_rad.cache[1]) == length(vlodf_rad.axes[2])
end
