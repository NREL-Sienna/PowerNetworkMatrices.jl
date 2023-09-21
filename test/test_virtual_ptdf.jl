@testset "Test Virtual PTDF matrices" begin
    sys = PSB.build_system(PSB.PSYTestSystems, "tamu_ACTIVSg2000_sys")
    ptdf_complete = PTDF(sys; linear_solver = "KLU")
    ptdf_virtual = VirtualPTDF(sys)

    for i in axes(ptdf_complete, 2)
        virtual = ptdf_virtual[i, :]
        for j in axes(ptdf_complete, 1)
            # check values using PTDFs axes
            @test isapprox(ptdf_complete[i, j], ptdf_virtual[i, j]; atol = 1e-10)
        end
    end
    # Check the cache is populated
    @test length(ptdf_virtual.cache) == length(ptdf_virtual.axes[1])
end

@testset "Test Virtual PTDF matrices with tolerance" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")
    ptdf_reference = deepcopy(S14_slackB1)
    ptdf_reference[abs.(ptdf_reference) .<= 1e-2] .= 0
    ptdf_virtual_with_tol = VirtualPTDF(sys; tol = 1e-2)
    for (n, i) in enumerate(axes(ptdf_virtual_with_tol, 1))
        # get the row
        @test isapprox(
            ptdf_virtual_with_tol[i, :], ptdf_reference[n, :], atol = 1e-3)
    end

    @test isapprox(
        sum(abs.(ptdf_reference[ptdf_virtual_with_tol.lookup[1]["Line12"], :])),
        sum(abs.(ptdf_virtual_with_tol["Line12", :])),
        atol = 1e-4)
end

@testset "Test Virtual PTDF matrices for 10 bus system with 2 reference buses" begin
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

@testset "Test Virtual PTDF cache" begin
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

@testset "Test Virtual PTDF with distributed slack" begin
    # get 5 bus system
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")
    bus_number = length(PNM.get_buses(sys))
    dist_slack = 1 / bus_number * ones(bus_number)
    # compute full PTDF
    ptdf = PTDF(sys; dist_slack = dist_slack)
    # compute each row of the virtual PTDF and compare values
    vptdf = VirtualPTDF(sys; dist_slack = dist_slack)
    for row in 1:size(ptdf.data, 2)
        # evaluate the column (just needs one element)
        vptdf[row, 1]
        @test isapprox(vptdf.cache[row], ptdf[row, :], atol = 1e-5)
    end
end

@testset "Test Virtual PTDF matrix with distributed bus and with 2 reference buses" begin
    # check if an error is correctly thrown

    # 2 reference bus system
    sys = PSB.build_system(PSISystems, "2Area 5 Bus System")
    buscount = length(PNM.get_buses(sys))
    dist_slack = 1 / buscount * ones(buscount)
    slack_array = dist_slack / sum(dist_slack)

    ptdf_1 = VirtualPTDF(sys; dist_slack = slack_array)
    @test_throws ErrorException ptdf_1[1, 1]

    # incorrect dist_slack array length
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    buscount = length(PNM.get_buses(sys5)) + 1
    dist_slack = 1 / buscount * ones(buscount)
    slack_array = dist_slack / sum(dist_slack)
    test_val2 = false
    ptdf_2 = VirtualPTDF(sys5; dist_slack = slack_array)
    @test_throws ErrorException ptdf_2[1, 1]
end

@testset "Test Virtual PTDF auxiliary functions" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")

    # test empty cache
    vptdf = VirtualPTDF(sys)
    @test isempty(vptdf) == true

    # test full cache
    vptdf[2, 1]
    @test isempty(vptdf) == false

    # check if error is correctly thrown
    @test_throws ErrorException vptdf[1, 1] = 1

    # get the rows and full PTDF matrix, test get_ptdf_data
    ptdf = PTDF(sys)
    for i in PNM.get_branch_ax(vptdf)
        for j in PNM.get_bus_ax(vptdf)
            vptdf[i, j]
        end
    end
    dict_ = Dict()
    for (n, i) in enumerate(PNM.get_branch_ax(vptdf))
        dict_[n] = vptdf.cache[n]
    end
    @test get_ptdf_data(vptdf) == dict_

    # test get axes values
    @test setdiff(PNM.get_branch_ax(vptdf), PSY.get_name.(PNM.get_ac_branches(sys))) ==
          String[]
    @test setdiff(PNM.get_bus_ax(vptdf), PSY.get_number.(PNM.get_buses(sys))) == String[]

    # test show
    test_value = false
    try
        show(vptdf)
        test_value = true
    catch err
        if err isa Exception
            test_value = false
        end
    end
    @test test_value
end
