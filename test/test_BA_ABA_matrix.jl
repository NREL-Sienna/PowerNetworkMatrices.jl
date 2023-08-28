@testset "Test A, BA, and ABA matrix creation" begin
    # test on 5 and 14 bus system
    for name in ["c_sys5", "c_sys14"]
        sys = PSB.build_system(PSB.PSITestSystems, name)
        # at first let's see if factorization flag works
        ABA_no_lu = ABA_Matrix(sys)
        @test isnothing(ABA_no_lu.K)
        # check the is_factorized function
        @test is_factorized(ABA_no_lu) == false
        # factorize the ABA matrix
        ABA_no_lu = factorize(ABA_no_lu)
        @test is_factorized(ABA_no_lu) == true
        # get the ABA matrix with the current method
        ABA_lu = ABA_Matrix(sys; factorize = true)
        # check the is_factorized function
        @test is_factorized(ABA_no_lu) == true
        # evaluate if the ABA matrix is correct
        A = IncidenceMatrix(sys)
        BA = BA_Matrix(sys)
        ABA_2 = transpose(A.data) * BA.data'
        @test isapprox(
            ABA_lu.data,
            ABA_2[setdiff(1:end, A.ref_bus_positions), setdiff(1:end, A.ref_bus_positions)],
            atol = 1e-8,
        )

        # evaluate if the LU factorization evaluates a correct PTDF matrix
        ptdf_1 = PTDF(sys)
        Ix = Matrix(1.0I, size(ABA_lu.data, 1), size(ABA_lu.data, 1))
        ABA_inv = ABA_lu.K \ Ix
        ptdf_2 = BA.data[setdiff(1:end, A.ref_bus_positions), :]' * ABA_inv
        @test isapprox(
            ptdf_1.data[setdiff(1:end, A.ref_bus_positions), :],
            ptdf_2',
            atol = 1e-8,
        )
    end
end

@testset "Test BA and ABA matrix indexing" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")

    # get the matrices
    ba = BA_Matrix(sys)
    aba = ABA_Matrix(sys)

    # check if indexing for the BA is correct (row and column indices)
    # ba matrix is stored as transposed
    for i in 1:size(ba, 2)
        for j in 1:size(ba, 1)
            @test isapprox(ba[i, j], ba.data[j, i])
        end
    end
    # check if indexing for the BA is correct (line names and bus numbers)
    lookup1 = PNM.get_lookup(ba)
    for i in axes(ba, 2)
        for j in axes(ba, 1)
            @test isapprox(ba[i, j], ba.data[lookup1[1][j], lookup1[2][i]])
        end
    end

    # check indexing for the ABA matrix
    lookup2 = aba.lookup[1]
    for i in axes(aba, 1)
        @test aba[i, :] == aba.data[lookup2[i], :]
    end

    # test if error is correctly thrown when ref bus is called
    rb = collect(ba.ref_bus_positions)[1]
    test_val = false
    try
        aba[rb, :]
    catch err
        if err isa ErrorException
            test_val = true
        else
            error("Expected an ErrorException but was not thrown.")
        end
    end
    @test test_val
end

@testset "Test show for A, BA and ABA matrix" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")

    a = IncidenceMatrix(sys)
    ba = BA_Matrix(sys)
    aba = ABA_Matrix(sys)

    for mat in [a, ba, aba]
        test_value = false
        try
            show(@eval a)
            test_value = true
        catch err
            if err isa Exception
                test_value = false
            end
        end
        @test test_value
    end
end
