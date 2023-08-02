@testset "Virtual PTDF matrices" begin
    sys = PSB.build_system(PSB.PSYTestSystems, "tamu_ACTIVSg2000_sys")
    ptdf_complete = PTDF(sys; linear_solver = "KLU")
    ptdf_virtual = VirtualPTDF(sys)

    for i in axes(ptdf_complete, 2)
        comp = ptdf_complete[i, :]
        virtual = ptdf_virtual[i, :]
        for j in axes(ptdf_complete, 1)
            # check values using PTDFs axes
            @test isapprox(ptdf_complete[i, j], ptdf_virtual[i, j]; atol = 1e-10)
        end
    end
    # Check the cache is populated
    @test length(ptdf_virtual.cache) == length(ptdf_virtual.axes[1])
end

@testset "Virtual PTDF matrices with tolerance" begin
    sys = PSB.build_system(PSB.PSYTestSystems, "tamu_ACTIVSg2000_sys")
    ptdf_virtual = VirtualPTDF(sys)
    ptdf_virtual_with_tol = VirtualPTDF(sys; tol = 1e-2)
    @test sum(abs.(ptdf_virtual["ODESSA 2 0  -1001-ODESSA 3 0  -1064-i_1", :])) >
          sum(abs.(ptdf_virtual_with_tol["ODESSA 2 0  -1001-ODESSA 3 0  -1064-i_1", :]))
end

@testset "Virtual PTDF matrices for 10 bus system with 2 reference buses" begin
    # get system
    sys = PSB.build_system(PSISystems, "2Area 5 Bus System")   # get the system composed by 2 5-bus ones connected by a DC line
    # get PTDF matrix with KLU as reference
    ptdf_complete = PTDF(sys; linear_solver = "KLU")
    # check VirtualPTDF rows with the ones from KLU
    ptdf_virtual = VirtualPTDF(sys)
    for i in axes(ptdf_complete, 2)
        comp = ptdf_complete[i, :]
        virtual = ptdf_virtual[i, :]
        # check values using PTDFs axes
        @test isapprox(comp, virtual; atol = 1e-10)
    end

    # check submatrices: since connected by a single bus, areas must have the same numbers
    branch_number = length(ptdf_complete.axes[2])
    bus_number = length(ptdf_complete.axes[1])
    ptdf_first_area = zeros(Int(branch_number / 2), Int(bus_number / 2))
    ptdf_second_area = zeros(Int(branch_number / 2), Int(bus_number / 2))

    for i in 1:Int(branch_number / 2)
        aa = ptdf_virtual.cache[i]
        bb = ptdf_virtual.cache[i + Int(branch_number / 2)]
        ptdf_first_area[i, :] .= aa[1:Int(bus_number / 2)]
        ptdf_second_area[i, :] .= bb[(Int(bus_number / 2) + 1):end]
    end
    @test isapprox(ptdf_first_area, ptdf_second_area, atol = 1e-6)
end

@testset "Test virtual PTDF cache" begin
    RTS = build_system(PSITestSystems, "test_RTS_GMLC_sys")
    line_names = get_name.(get_components(Line, RTS))
    persist_lines = line_names[1:10]

    vptdf = VirtualPTDF(RTS; max_cache_size = 1, persistent_lines = persist_lines)

    for l in line_names
        @test size(vptdf[l, :]) == (73,)
    end

    for l in persist_lines
        @test vptdf.lookup[1][l] ∈ keys(vptdf.cache.temp_cache)
    end
end

@testset "Test virtual PTDF with distributed slack" begin
    # get 5 bus system
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    # compute full PTDF
    ptdf = PTDF(sys; dist_slack = slack_array)
    # compute each row of the virtual PTDF and compare values
    vptdf = VirtualPTDF(sys; dist_slack = slack_array)
    for row in 1:size(ptdf.data, 2)
        # evaluate the column (just needs one element)
        vptdf[row, 1]
        @assert isapprox(vptdf.cache[row], ptdf[row, :], atol = 1e-5)
    end
end
