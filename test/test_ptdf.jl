@testset "Test PTDF matrices, w/ and w/o tolerance" begin
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    for solver in ["KLU", "Dense", "MKLPardiso"]
        for approach in ["standard", "separate_matrices"]
            buses_5 = nodes5()
            branches_5 = branches5(buses_5)
            if approach == "standard"
                P5 = PTDF(branches_5, buses_5; linear_solver = solver)
                P5_bis = PTDF(sys5; linear_solver = solver)
                P5_sparse = PTDF(branches_5, buses_5; linear_solver = solver, tol = 1e-3)
            elseif approach == "separate_matrices"
                if solver == "Dense"
                    continue
                else
                    A = IncidenceMatrix(sys5)
                    BA = BA_Matrix(sys5)
                    P5 = PTDF(A, BA; linear_solver = solver)
                    P5_sparse = PTDF(A, BA; linear_solver = solver, tol = 1e-3)
                end
            end

            # test method with branches and buses
            @test isapprox(maximum(P5.data), maximum(S5_slackB4), atol = 1e-3)
            @test isapprox(P5[branches_5[1], buses_5[1]], 0.1939166051164976)
            # test method with whole system
            if approach == "standard"
                @test isapprox(maximum(P5_bis.data), maximum(S5_slackB4), atol = 1e-3)
                @test isapprox(P5_bis[branches_5[1], buses_5[1]], 0.1939166051164976)
            end

            # additional test with other set of branches and buses
            buses_14 = nodes14()
            branches_14 = branches14(buses_14)
            P14 = PTDF(branches_14, buses_14)
            @test isapprox(maximum(P14.data), maximum(S14_slackB1), atol = 1e-3)

            # check getindex (line name and bus number)
            P5NS = PTDF([branches_5[b] for b in Br5NS_ids], [buses_5[b] for b in Bu5NS_ids])
            for br in Br5NS_ids, b in Bu5NS_ids
                @test isapprox(
                    getindex(P5NS, string(br), b),
                    S5_slackB4[br, b],
                    atol = 1e-3,
                )
            end

            # check getindex (iterated accoring to rows and columns)
            P5NS_mod = zeros(size(S5_slackB4))
            for br in Br5NS_ids, b in Bu5NS_ids
                # get the correct indices
                row = P5NS.lookup[2][string(br)]
                col = P5NS.lookup[1][b]
                P5NS_mod[row, col] = S5_slackB4[br, b]
                @test isapprox(
                    getindex(P5NS, row, col),
                    S5_slackB4[br, b],
                    atol = 1e-3,
                )
            end

            # adjust the matrix to have same indexing
            @test isapprox(get_ptdf_data(P5NS), P5NS_mod, atol = 1e-3)

            PRTS = PTDF(RTS; linear_solver = solver)
            PRTS_mod = zeros(size(SRTS_GMLC))
            PRTS_sparse = PTDF(RTS; linear_solver = solver, tol = 1e-3)
            bnums = sort([PSY.get_number(b) for b in PSY.get_components(Bus, RTS)])
            for (ibr, br) in enumerate(RTS_branchnames), (ib, b) in enumerate(bnums)
                row = PRTS.lookup[2][string(br)]
                col = PRTS.lookup[1][b]
                PRTS_mod[row, col] = SRTS_GMLC[ibr, ib]
                @test isapprox(getindex(PRTS, br, b), SRTS_GMLC[ibr, ib], atol = 1e-3)
            end

            # manually sparsify the matrix
            PRTS.data[abs.(PRTS.data) .< 1e-3] .= 0
            @test sum(abs.(get_ptdf_data(PRTS) - get_ptdf_data(PRTS_sparse)) .> 1e-3) == 0
        end
    end

    # check axes values
    P5 = PTDF(sys5)
    @test setdiff(PNM.get_branch_ax(P5), PSY.get_name.(PNM.get_ac_branches(sys5))) ==
          String[]
    @test setdiff(PNM.get_bus_ax(P5), PSY.get_number.(PNM.get_buses(sys5))) == String[]

    # auxiliary function
    PRTS_sparse = PTDF(RTS; tol = 1e-3)
    @test PNM.get_tol(PRTS_sparse).x == Base.RefValue(1e-3).x
end

@testset "Test PTDF matrices for 10 bus system with 2 reference buses" begin
    # get system
    sys = PSB.build_system(PSISystems, "2Area 5 Bus System")   # get the system composed by 2 5-bus ones connected by a DC line
    ptdf_complete_klu = PTDF(sys; linear_solver = "KLU")
    ptdf_complete_dense = PTDF(sys; linear_solver = "Dense")

    @test sum(ptdf_complete_klu.data - ptdf_complete_dense.data) < 1e-9
    @test isapprox(ptdf_complete_klu.data, ptdf_complete_dense.data, atol = 1e-6)

    # check submatrices: siunce connected by a single bus, areas must have the same numbers
    branch_number = length(ptdf_complete_klu.axes[1])
    bus_number = length(ptdf_complete_klu.axes[2])
    ptdf_first_area =
        ptdf_complete_klu.data[1:Int(branch_number / 2), 1:Int(bus_number / 2)]
    ptdf_second_area = ptdf_complete_klu.data[
        (Int(branch_number / 2) + 1):end,
        (Int(bus_number / 2) + 1):end,
    ]
    @test isapprox(ptdf_first_area, ptdf_second_area, atol = 1e-6)
end

@testset "Test serialization of PTDF matrices to HDF5" begin
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    buses_5 = nodes5()
    branches_5 = branches5(buses_5)
    P5 = PTDF(branches_5, buses_5; linear_solver = "KLU")
    P5_sparse = PTDF(branches_5, buses_5; linear_solver = "KLU", tol = 1e-3)
    for ptdf in (P5, P5_sparse)
        for compress in (true, false)
            path = mktempdir()
            filename = joinpath(path, "ptdf.h5")
            @test !isfile(filename)
            to_hdf5(ptdf, filename; compress = compress)
            @test isfile(filename)
            ptdf2 = PTDF(filename)
            @test ptdf == ptdf2
        end
    end
end

@testset "Test System with isolated buses" begin
    sys_1 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    PSY.add_component!(
        sys_1,
        PSY.ACBus(;
            number = 6,
            name = "isolated_node_1",
            bustype = PSY.ACBusTypes.ISOLATED,
            angle = 0.0,
            magnitude = 1.1,
            voltage_limits = (min = 0.9, max = 1.1),
            base_voltage = 230.0,
        ),
    )
    PSY.add_component!(
        sys_1,
        PSY.ACBus(;
            number = 7,
            name = "isolated_node_2",
            bustype = PSY.ACBusTypes.ISOLATED,
            angle = 0.0,
            magnitude = 1.1,
            voltage_limits = (min = 0.9, max = 1.1),
            base_voltage = 230.0,
        ),
    )
    ptdf_1 = PTDF(sys_1)
    # Test that the isolated buses are not included in the PTDF matrix
    @test length(axes(ptdf_1)[1]) == 5

    sys_2 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    branches_2 = PNM.get_ac_branches(sys_2)

    PSY.add_component!(
        sys_2,
        PSY.ACBus(;
            number = 6,
            name = "isolated_node_1",
            bustype = PSY.ACBusTypes.ISOLATED,
            angle = 0.0,
            magnitude = 1.1,
            voltage_limits = (min = 0.9, max = 1.1),
            base_voltage = 230.0,
        ),
    )

    add_component!(
        sys_2,
        PSY.Line(;
            name = "7",
            available = branches_2[2].available,
            active_power_flow = branches_2[2].active_power_flow,
            reactive_power_flow = branches_2[2].reactive_power_flow,
            arc = PSY.Arc(;
                from = PSY.get_component(PSY.ACBus, sys_2, "nodeA"),
                to = PSY.get_component(PSY.ACBus, sys_2, "isolated_node_1"),
            ),
            r = branches_2[2].r,
            x = branches_2[2].x,
            b = branches_2[2].b,
            rating = get_rating(branches_2[2]),
            angle_limits = get_angle_limits(branches_2[2]),
        ),
    )

    # Test Throw error when isolated buses are connected to an available branch
    @test_throws IS.ConflictingInputsError ptdf_2 = PTDF(sys_2)
end

@testset "Test PTDF matrices with distributed slack" begin
    """
    CAUTION: this test just test that all the matrices are the same, but there
    are no reference values.
    """

    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")

    buscount = length(PNM.get_buses(sys5))

    dist_slack = 1 / buscount * ones(buscount)
    slack_array = dist_slack / sum(dist_slack)

    P5_1 = PTDF(sys5; dist_slack = slack_array, linear_solver = "KLU")
    P5_2 = PTDF(sys5; dist_slack = slack_array, linear_solver = "Dense")
    P5_3 = PTDF(sys5; dist_slack = slack_array, linear_solver = "MKLPardiso")

    @test isapprox(P5_1.data, P5_2.data, atol = 1e-5)
    @test isapprox(P5_1.data, P5_3.data, atol = 1e-5)
    @test isapprox(P5_2.data, P5_3.data, atol = 1e-5)
end

@testset "Test PTDF matrix with distributed bus and with 2 reference buses" begin

    # check if an error is correctly thrown

    # 2 reference bus system
    sys = PSB.build_system(PSISystems, "2Area 5 Bus System")
    buscount = length(PNM.get_buses(sys))
    dist_slack = 1 / buscount * ones(buscount)
    slack_array = dist_slack / sum(dist_slack)

    @test_throws ErrorException ptdf_1 = PTDF(sys; dist_slack = slack_array)

    # incorrect dist_slack arrya length
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    buscount = length(PNM.get_buses(sys5)) + 1
    dist_slack = 1 / buscount * ones(buscount)
    slack_array = dist_slack / sum(dist_slack)

    @test_throws ErrorException ptdf_2 = PTDF(sys5; dist_slack = slack_array)
end
