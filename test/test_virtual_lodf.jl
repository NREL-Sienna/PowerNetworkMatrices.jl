@testset "Test Virtual LODF matrices" begin
    sys = PSB.build_system(PSB.PSYTestSystems, "tamu_ACTIVSg2000_sys")
    vlodf = VirtualLODF(sys)
    LODF_ref = LODF(sys)

    # data
    for (idx, name) in enumerate(vlodf.axes[1])
        # compare
        @test isapprox(vlodf[name, :], LODF_ref[idx, :], atol = 1e-6)
    end

    # check dicionary with rows
    data_dict = get_lodf_data(vlodf)
    for i in 1:length(vlodf.axes[1])
        # compare
        @test isapprox(data_dict[i], LODF_ref[i, :], atol = 1e-6)
    end

    # check the getindex function works properly
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    vlodf5 = VirtualLODF(sys5)
    for i in axes(Lodf_5, 1)
        for j in axes(Lodf_5, 2)
            # get the data
            @test isapprox(vlodf5[i, j], Lodf_5[i, j], atol = 1e-3)
        end
    end
end

@testset "Test Virtual LODF functions" begin
    # get system
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    buses_5 = nodes5()
    branches_5 = branches5(buses_5)
    vlodf = VirtualLODF(sys5)
    # properties
    @test size(vlodf) == (length(branches_5), length(branches_5))
    @test length(PNM.get_branch_ax(vlodf)) == length(branches_5)
end

@testset "Test Virtual LODF matrices with tolerance" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")
    lodf_reference = deepcopy(Lodf_14)
    lodf_reference[abs.(lodf_reference) .<= 1e-2] .= 0
    lodf_virtual_with_tol = VirtualLODF(sys; tol = 1e-2)
    for (n, i) in enumerate(axes(lodf_virtual_with_tol, 1))
        # get the row
        @test isapprox(
            lodf_virtual_with_tol[i, :], lodf_reference[n, :], atol = 1e-3)
    end

    @test isapprox(
        sum(abs.(lodf_reference[lodf_virtual_with_tol.lookup[1]["Line12"], :])),
        sum(abs.(lodf_virtual_with_tol["Line12", :])),
        atol = 1e-5)
end

@testset "Test Virtual LODF matrices for 10 bus system with 2 reference buses" begin
    # get system
    sys = PSB.build_system(PSISystems, "2Area 5 Bus System")   # get the system composed by 2 5-bus ones connected by a DC line
    # get PTDF matrix with KLU as reference
    lodf_complete = LODF(sys; linear_solver = "KLU")
    # check VirtualPTDF rows with the ones from KLU
    lodf_virtual = VirtualLODF(sys)
    for i in axes(lodf_complete, 2)
        comp = lodf_complete[i, :]
        virtual = lodf_virtual[i, :]
        # check values using PTDFs axes
        @test isapprox(comp, virtual; atol = 1e-10)
    end

    # check submatrices: since connected by a single bus, areas must have the same numbers
    branch_number = length(lodf_complete.axes[2])
    bus_number = length(lodf_complete.axes[1])
    ptdf_first_area = zeros(Int(branch_number / 2), Int(bus_number / 2))
    ptdf_second_area = zeros(Int(branch_number / 2), Int(bus_number / 2))

    for i in 1:Int(branch_number / 2)
        aa = lodf_virtual.cache[i]
        bb = lodf_virtual.cache[i + Int(branch_number / 2)]
        ptdf_first_area[i, :] .= aa[1:Int(bus_number / 2)]
        ptdf_second_area[i, :] .= bb[(Int(bus_number / 2) + 1):end]
    end
    @test isapprox(ptdf_first_area, ptdf_second_area, atol = 1e-6)
end

@testset "Test virtual LODF cache" begin
    RTS = build_system(PSITestSystems, "test_RTS_GMLC_sys")
    line_names = get_name.(PNM.get_ac_branches(RTS))
    persist_lines = line_names[1:10]

    vlodf = VirtualLODF(RTS; max_cache_size = 1, persistent_lines = persist_lines)

    for l in line_names
        @test size(vlodf[l, :]) == (120,)
    end

    for l in persist_lines
        @test vlodf.lookup[1][l] ∈ keys(vlodf.cache.temp_cache)
    end
end
