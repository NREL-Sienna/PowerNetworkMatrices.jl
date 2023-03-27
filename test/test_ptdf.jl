@testset "PTDF matrices" begin
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    for solver in ["KLU", "Dense", "MKLPardiso"]
        for approach in ["standard", "separate_matrices"]
            nodes_5 = nodes5()
            branches_5 = branches5(nodes_5)
            if approach == "standard"
                P5 = PTDF(branches_5, nodes_5; linear_solver = solver)
            elseif approach == "separate_matrices"
                if solver == "Dense"
                    continue
                else
                    A = IncidenceMatrix(sys5)
                    BA = BA_Matrix(sys5)
                    P5 = PTDF(A, BA; linear_solver = solver)
                end
            end
            @test isapprox(maximum(P5.data), maximum(S5_slackB4), atol = 1e-3)
            @test isapprox(P5[branches_5[1], nodes_5[1]], 0.1939166051164976)

            nodes_14 = nodes14()
            branches_14 = branches14(nodes_14)
            P14 = PTDF(branches_14, nodes_14)
            @test isapprox(maximum(P14.data), maximum(S14_slackB1), atol = 1e-3)

            P5NS = PTDF([branches_5[b] for b in Br5NS_ids], [nodes_5[b] for b in Bu5NS_ids])
            for br in Br5NS_ids, b in Bu5NS_ids
                @test isapprox(
                    getindex(P5NS, string(br), b),
                    S5_slackB4[br, b],
                    atol = 1e-3,
                )
            end

            PRTS = PTDF(RTS)
            bnums = sort([PSY.get_number(b) for b in PSY.get_components(Bus, RTS)])
            for (ibr, br) in enumerate(RTS_branchnames), (ib, b) in enumerate(bnums)
                @test isapprox(getindex(PRTS, br, b), SRTS_GMLC[ibr, ib], atol = 1e-3)
            end
        end
    end
end

@testset "PTDF matrices for 10 bus system with 2 reference buses" begin
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
