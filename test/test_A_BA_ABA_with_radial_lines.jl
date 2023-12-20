# TODO: missing tests for BA matrix with radial lines
# TODO: missing tests for ABA matrix with radial lines (creating ABA from A and BA with radial lines)
@testset "Test A matrix with radial lines" begin
    # A matrix evaluated with radial lines must be the same as the ogirinal A
    # matrix without the rows and columns related to radiale branches and leaf nodes
    for name in ["c_sys14", "test_RTS_GMLC_sys"]
        # load the system
        sys = PSB.build_system(PSB.PSITestSystems, name);
        # get the incidence matrix
        A = IncidenceMatrix(sys);
        # ... and with radial lines
        A_rad = IncidenceMatrix(sys; reduce_radial_branches=true)
        # get inidices for the leaf nodes
        rb = A_rad.reduce_radial_branches
        bus_numbers = []
        for i in keys(rb.bus_reduction_map)
            append!(bus_numbers, collect(rb.bus_reduction_map[i]))
        end
        bus_idx = setdiff(1:size(A.data, 2), [A.lookup[2][i] for i in bus_numbers])
        # ... and radial branches
        br_idx = setdiff(1:size(A.data, 1), [A.lookup[1][i] for i in rb.radial_branches])
        # now extract A matrix anc compare
        @test all(isapprox.(A.data[br_idx, bus_idx], A_rad.data))
    end
end

@testset "Test BA matrix with radial lines" begin
    for name in ["c_sys14", "test_RTS_GMLC_sys"]
        # load the system
        sys = PSB.build_system(PSB.PSITestSystems, name);
        # get the incidence matrix
        BA = BA_Matrix(sys);
        # ... and with radial lines
        BA_rad = IncidenceMatrix(sys; reduce_radial_branches=true);
        # get inidices for the leaf nodes
        rb = BA_rad.reduce_radial_branchesa;
        bus_numbers = [];
        for i in keys(rb.bus_reduction_map)
            append!(bus_numbers, collect(rb.bus_reduction_map[i]))
        end
        bus_idx = setdiff(1:size(A.data, 2), [A.lookup[2][i] for i in bus_numbers])
        # ... and radial branches
        br_idx = setdiff(1:size(A.data, 1), [A.lookup[1][i] for i in rb.radial_branches])
        # now extract A matrix anc compare
        @test all(isapprox.(BA.data[bus_idx, br_idx], BA_rad.data))
    end
end

@testset "Test ABA matrix with radial lines" begin
    for name in ["c_sys14", "test_RTS_GMLC_sys"]
        # load the system
        sys = PSB.build_system(PSB.PSITestSystems, name)
        # at first, get the radial branches
        rb = RadialBranches(IncidenceMatrix(sys))
        # get the ABA matrix
        ABA_lu = ABA_Matrix(sys)
        # evaluate if the ABA matrix is correct
        A = IncidenceMatrix(sys)
        BA = BA_Matrix(sys)
        ABA_2 = transpose(A.data) * BA.data'
        @test isapprox(
            ABA_lu.data,
            ABA_2[setdiff(1:end, A.ref_bus_positions), setdiff(1:end, A.ref_bus_positions)],
            atol = 1e-8,
        )
    end
end