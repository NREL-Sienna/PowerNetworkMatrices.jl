@testset "Test ABA matrix" begin
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

const PNM = PowerNetworkMatrices
branches = PNM.get_ac_branches(sys)
nodes = PNM.get_buses(sys)

line_ax = [PSY.get_name(branch) for branch in branches]
bus_ax = [PSY.get_number(bus) for bus in nodes]
axes = (line_ax, bus_ax)
M, bus_ax_ref = PNM.calculate_adjacency(branches, nodes)
ref_bus_positions = PNM.find_slack_positions(nodes)
look_up = (PNM.make_ax_ref(line_ax), bus_ax_ref)
dist_slack = Float64[]

A, ref_bus_positions = PNM.calculate_A_matrix(branches, nodes, ref_bus_positions)
BA = PNM.calculate_BA_matrix(branches, look_up[2])

linecount = size(BA, 2)
buscount = size(BA, 1)

using KLU
using LinearAlgebra

ABA = PNM.calculate_ABA_matrix(A, BA, ref_bus_positions)
K = klu(ABA)
# inizialize matrices for evaluation
PTDFm_t = zeros(buscount, linecount)
valid_ix = setdiff(1:size(PTDFm_t, 1), ref_bus_positions)
copyto!(PTDFm_t, BA)
KLU.solve!(K, PTDFm_t[valid_ix, :])
