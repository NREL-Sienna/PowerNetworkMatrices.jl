@testset "Test LODF matrices" begin
    """
    NOTE: LODF is transposed
    """

    # get 5 bus system
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    buses_5 = nodes5()
    branches_5 = branches5(buses_5)

    # get LODF with buses and branches
    L5 = LODF(branches_5, buses_5)
    @test isapprox(maximum(L5.data), maximum(Lodf_5), atol = 1e-3)
    @test isapprox(L5[branches_5[1], branches_5[2]], 0.3447946513849093)

    # get LODF with second method
    a = IncidenceMatrix(sys5)
    ptdf = PTDF(sys5)
    lodf_t_3 = PNM._calculate_LODF_matrix_KLU2(a.data, ptdf.data)
    @test isapprox(lodf_t_3, L5.data, atol = 1e-5)

    # get LODF with system
    L5NS = LODF(sys5)
    @test getindex(L5NS, "5", "6") - -0.3071 <= 1e-4
    total_error = abs.(L5NS.data' .- Lodf_5)
    @test isapprox(sum(total_error), 0.0, atol = 1e-3)

    # A and PTDF case
    A = IncidenceMatrix(sys5)
    P5 = PTDF(sys5)
    L5NS_from_ptdf = LODF(A, P5)
    @test getindex(L5NS_from_ptdf, "5", "6") - -0.3071 <= 1e-4
    total_error = abs.(L5NS_from_ptdf.data' .- Lodf_5)
    @test isapprox(sum(total_error), 0.0, atol = 1e-3)

    # A, ABA, and BA case
    ABA = ABA_Matrix(sys5; factorize = true)
    BA = BA_Matrix(sys5)
    L5NS_from_ba_aba = LODF(A, ABA, BA)
    @test getindex(L5NS_from_ba_aba, "5", "6") - -0.3071 <= 1e-4
    total_error = abs.(L5NS_from_ba_aba.data' .- Lodf_5)
    @test isapprox(sum(total_error), 0.0, atol = 1e-3)

    for i in axes(Lodf_5, 1)
        @test isapprox(L5[i, :], Lodf_5[i, :], atol = 1e-3)
        @test isapprox(L5NS[i, :], Lodf_5[i, :], atol = 1e-3)
        @test isapprox(L5NS_from_ptdf[i, :], Lodf_5[i, :], atol = 1e-3)
        @test isapprox(L5NS_from_ba_aba[i, :], Lodf_5[i, :], atol = 1e-3)
    end

    # test if error is thrown in case `tol` is defined in PTDF
    P5 = PTDF(sys5; tol = 1e-3)
    test_value = false
    try
        L5NS_from_ptdf = LODF(A, P5)
    catch err
        if err isa ErrorException
            test_value = true
        end
    end
    @test test_value

    # test if error is thrown in case `MKLPardiso` is chosen as a linera solver
    test_value = false
    try
        lodf = LODF(sys5; linear_solver = "MKLPardiso")
    catch err
        if err isa ErrorException
            test_value = true
        end
    end
    @test test_value

    # get 14 bus system
    buses_14 = nodes14()
    branches_14 = branches14(buses_14)
    L14 = LODF(branches_14, buses_14)
    @test isapprox(maximum(L14.data), maximum(Lodf_14), atol = 1e-3)
end

@testset "Test sparse LODF matrix" begin
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")

    # basic case: system
    L5NS_1 = LODF(sys5; tol = 0.4)

    # A and PTDF case
    a_matrix = IncidenceMatrix(sys5)
    ptdf_matrix = PTDF(sys5)
    L5NS_3 = LODF(a_matrix, ptdf_matrix; tol = 0.4)

    # A, ABA, and BA case
    aba_matrix = ABA_Matrix(sys5; factorize = true)
    ba_matrix = BA_Matrix(sys5)
    L5NS_4 = LODF(a_matrix, aba_matrix, ba_matrix; tol = 0.4)

    # reference value
    ref_sparse_Lodf_5 = deepcopy(Lodf_5')
    ref_sparse_Lodf_5[abs.(ref_sparse_Lodf_5) .< PNM.get_tol(L5NS_1)] .= 0

    # tests
    @test isapprox(Matrix(L5NS_1.data), ref_sparse_Lodf_5, atol = 1e-3)
    @test isapprox(Matrix(L5NS_1.data), L5NS_3.data, atol = 1e-3)
    @test isapprox(Matrix(L5NS_1.data), L5NS_4.data, atol = 1e-3)
end

@testset "Test LODF getindex and get_lodf_data" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")

    # get the LODF matrix
    lodf_ = LODF(sys)

    # test get_lodf_data
    lodf_t = lodf_.data'
    @test isapprox(lodf_t, get_lodf_data(lodf_))

    # test getindex
    for name1 in PNM.get_branch_ax(lodf_)       # selected line
        for name2 in PNM.get_branch_ax(lodf_)   # outage line
            i = lodf_.lookup[1][name1]
            j = lodf_.lookup[2][name2]
            element_1 = lodf_[name1, name2]
            element_2 = lodf_[i, j]
            element_3 = lodf_t[i, j]
            @test isapprox(element_1, element_2, atol = 1e-5)
            @test isapprox(element_1, element_3, atol = 1e-5)
            @test isapprox(element_2, element_3, atol = 1e-5)
        end
    end
end