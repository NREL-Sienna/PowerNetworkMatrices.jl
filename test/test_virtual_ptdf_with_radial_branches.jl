@testset "Test VirtualPTDF with radial lines" begin
    # get the system
    sys = PSB.build_system(PSB.PSITestSystems, "test_RTS_GMLC_sys")
    # get the PTDF matrix for reference
    ptdf_rad = PTDF(sys; reduce_radial_branches = true)
    vptdf_rad = VirtualPTDF(sys; reduce_radial_branches = true)

    for i in axes(ptdf_rad, 2)
        virtual = vptdf_rad[i, :]
        for j in axes(ptdf_rad, 1)
            # check values using PTDFs axes
            @test isapprox(ptdf_rad[i, j], vptdf_rad[i, j]; atol = 1e-10)
        end
    end
    # Check the cache is populated
    @test length(vptdf_rad.cache) == length(vptdf_rad.axes[1])
end

@testset "Test VirtualPTDF with radial lines and distributed slack" begin
    # check if VirtualPTDF have same values as PTDF row-wise
    # sys = PSB.build_system(PSB.PSITestSystems, "test_RTS_GMLC_sys")
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")   # ! remove
    buscount = length(PNM.get_buses(sys))
    dist_slack = 1 / buscount * ones(buscount)
    slack_array = dist_slack / sum(dist_slack)
    ptdf_rad = PTDF(sys; reduce_radial_branches = true, dist_slack = slack_array)
    vptdf_rad = VirtualPTDF(sys; reduce_radial_branches = true, dist_slack = slack_array)
    for i in axes(ptdf_rad, 2)
        virtual = vptdf_rad[i, :]
        for j in axes(ptdf_rad, 1)
            # check values using PTDFs axes
            @test isapprox(ptdf_rad[i, j], vptdf_rad[i, j]; atol = 1e-10)
        end
    end
end
