# if it fails, we don't want the terminal to be flooded with errors, therefore failfast=true
@testset failfast = true "Test Virtual PTDF matrices" begin
    sys =
        @test_logs (:error, r"no active generators found at bus") match_mode = :any PSB.build_system(
            PSB.PSYTestSystems,
            "tamu_ACTIVSg2000_sys",
        )
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

#TODO - need to define a mapping from the data 
@testset "Test Virtual PTDF matrices with tolerance" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")
    ix_to_arc_map = Dict()
    ix_to_bus_map = Dict()
    for (ix, br) in enumerate(get_components(ACBranch, sys))
        ix_to_arc_map[ix] = PNM.get_arc_tuple(br)
    end
    for (ix, b) in enumerate(get_components(ACBus, sys))
        ix_to_bus_map[ix] = get_number(b)
    end

    ptdf_reference = deepcopy(S14_slackB1)
    ptdf_virtual_with_tol = VirtualPTDF(sys; tol = 1e-2) #
    ptdf_with_tol = PTDF(sys; tol = 1e-2)
    ptdf_with_tol_rearranged = zeros(size(ptdf_reference)[2], size(ptdf_reference)[1])
    ptdf_virtual_with_tol_rearranged = zeros(size(ptdf_reference))
    for ix in 1:size(ptdf_reference)[1], jx in 1:size(ptdf_reference)[2]
        v_value = ptdf_with_tol[ix_to_arc_map[ix], ix_to_bus_map[jx]]
        @show ix, jx, v_value
        ptdf_with_tol_rearranged[jx, ix] = v_value
        #L5NS_3_rearranged[ix,jx]= L5NS_3[ix_to_tuple_map[jx],ix_to_tuple_map[ix] ]
        #L5NS_4_rearranged[ix,jx]= L5NS_4[ix_to_tuple_map[jx],ix_to_tuple_map[ix] ]
    end
    ptdf_with_tol_rearranged
    ptdf_virtual_with_tol_rearranged
    #ptdf_reference[abs.(ptdf_reference) .<= 1e-2] .= 0
    ptdf = PTDF(sys)
    for (n, i) in enumerate(axes(ptdf_virtual_with_tol, 1))
        # get the row
        @test isapprox(
            ptdf[ix_to_arc_map[n], :], ptdf_reference[n, :], atol = 1e-3)
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

    # check submatrices: system has identical areas connected by a single hvdc, areas must have the same numbers for corresponding buses and arcs
    corresponding_buses = [(1, 6), (2, 7), (3, 8), (4, 9), (5, 10)]
    corresponding_arcs = [
        ((1, 2), (6, 7)),
        ((4, 5), (9, 10)),
        ((1, 5), (6, 10)),
        ((1, 4), (6, 9)),
        ((3, 4), (8, 9)),
        ((2, 3), (7, 8)),
    ]
    for bus_pair in corresponding_buses, arc_pair in corresponding_arcs
        @test isapprox(
            ptdf_virtual[arc_pair[1], bus_pair[1]],
            ptdf_virtual[arc_pair[2], bus_pair[2]],
            atol = 1e-6,
        )
    end
end

#TODO - the ptdf of the RTS is missing an arc? 
@testset "Test Virtual PTDF cache" begin
    RTS = build_system(PSITestSystems, "test_RTS_GMLC_sys")
    arc_tuples = [PNM.get_arc_tuple(arc) for arc in get_components(Arc, RTS)]
    persist_arcs = arc_tuples[1:10]

    vptdf = VirtualPTDF(RTS; max_cache_size = 1, persistent_arcs = persist_arcs)
    for l in arc_tuples
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
    @test isempty(vptdf)

    # test full cache
    vptdf[2, 1]
    @test !isempty(vptdf)

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
    arc_tuples = [PNM.get_arc_tuple(arc) for arc in get_components(Arc, sys)]
    @test setdiff(PNM.get_branch_ax(vptdf), arc_tuples) == []
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
