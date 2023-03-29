@testset "Test ABA matrix" begin
    # test on 5 and 14 bus system
    for name in ["sys5", "sys14"]
        if name == "sys5"
            sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
        else
            sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")
        end
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
        ABA_2 = A.data[:, setdiff(1:end, A.ref_bus_positions)]' * BA.data
        @test isapprox(ABA_lu.data, ABA_2, atol = 1e-8)

        # evaluate if the LU factorization evaluates a correct PTDF matrix
        ptdf_1 = PTDF(sys)
        Ix = Matrix(1.0I, size(ABA_lu.data, 1), size(ABA_lu.data, 1))
        ABA_inv = ABA_lu.K \ Ix
        ptdf_2 = BA.data * ABA_inv
        @test isapprox(
            ptdf_1.data[:, setdiff(1:end, A.ref_bus_positions)],
            ptdf_2,
            atol = 1e-8,
        )
    end
end