@testset "LODF matrices" begin
    """
    NOTE: LODF is transposed
    """

    # get system
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    buses_5 = nodes5()
    branches_5 = branches5(buses_5)

    # get LODF
    L5 = LODF(branches_5, buses_5)
    @test isapprox(maximum(L5.data), maximum(Lodf_5), atol = 1e-3)
    @test isapprox(L5[branches_5[2], branches_5[1]], 0.3447946513849093)

    # get LODF with second method
    a = IncidenceMatrix(sys5)
    ptdf = PTDF(sys5)
    lodf_t_3 = PNM._calculate_LODF_matrix_KLU2(a.data, ptdf.data)
    @test isapprox(lodf_t_3, L5.data, atol = 1e-5)

    buses_14 = nodes14()
    branches_14 = branches14(buses_14)
    L14 = LODF(branches_14, buses_14)
    @test isapprox(maximum(L14.data), maximum(Lodf_14), atol = 1e-3)

    L5NS = LODF(sys5)
    @test getindex(L5NS, "6", "5") - -0.3071 <= 1e-4
    total_error = abs.(L5NS.data' .- Lodf_5)
    @test isapprox(sum(total_error), 0.0, atol = 1e-3)

    # A and PTDF case
    A = IncidenceMatrix(sys5)
    P5 = PTDF(sys5)
    L5NS_from_ptdf = LODF(A, P5)
    @test getindex(L5NS_from_ptdf, "6", "5") - -0.3071 <= 1e-4
    total_error = abs.(L5NS_from_ptdf.data' .- Lodf_5)
    @test isapprox(sum(total_error), 0.0, atol = 1e-3)

    # A, ABA, and BA case
    ABA = ABA_Matrix(sys5; factorize = true)
    BA = BA_Matrix(sys5)
    L5NS_from_ba_aba = LODF(A, ABA, BA)
    @test getindex(L5NS_from_ba_aba, "6", "5") - -0.3071 <= 1e-4
    total_error = abs.(L5NS_from_ba_aba.data' .- Lodf_5)
    @test isapprox(sum(total_error), 0.0, atol = 1e-3)
end

@testset "Sparse LODF matrix" begin
    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")

    # basic case: system
    L5NS_1 = LODF(sys5; tol = 0.4)

    # basic case: system but sparsification is done afterwards
    L5NS_2 = LODF(sys5)
    drop_small_entries!(L5NS_2, 0.4)

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
    ref_sparse_Lodf_5[abs.(ref_sparse_Lodf_5) .< 0.4] .= 0

    # tests
    @test isapprox(Matrix(L5NS_1.data), ref_sparse_Lodf_5, atol = 1e-3)
    @test isapprox(Matrix(L5NS_1.data), L5NS_2.data, atol = 1e-3)
    @test isapprox(Matrix(L5NS_1.data), L5NS_3.data, atol = 1e-3)
    @test isapprox(Matrix(L5NS_1.data), L5NS_4.data, atol = 1e-3)
end

@testset "LODF matrices with distributed slack" begin
    """
    CAUTION: this test just test that all the matrices are the same, but there 
    are no reference values.
    """

    sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")

    buscount = length(PNM.get_buses(sys5))

    dist_slack = 1 / buscount * ones(buscount)
    slack_array = dist_slack / sum(dist_slack)

    P5_1 = PTDF(sys5; dist_slack = slack_array, linear_solver = "KLU")

    # method with A and PTDF, PTDF already evaluated with distributed slack bus
    A = IncidenceMatrix(sys5)
    LODF_ref = LODF(A, P5_1)

    # base case with system data
    LODF_1 = LODF(sys5; dist_slack = slack_array)

    # tests
    @test isapprox(LODF_ref.data, LODF_1.data, atol = 1e-5)
end
