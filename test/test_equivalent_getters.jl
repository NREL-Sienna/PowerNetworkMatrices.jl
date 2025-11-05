@testset "Equivalent getters for BranchesParallel" begin
    # Create test system with parallel branches
    sys = PSB.build_system(PSB.PSITestSystems, "case10_radial_series_reductions")

    # Create two test lines with known parameters for parallel configuration
    bus1 = first(PSY.get_components(PSY.ACBus, sys))
    bus2 = collect(PSY.get_components(PSY.ACBus, sys))[2]

    # Create test branches with specific values
    line1 = PSY.Line(;
        name = "test_line_1",
        available = true,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = PSY.Arc(; from = bus1, to = bus2),
        r = 0.1,  # resistance
        x = 0.2,  # reactance
        b = (from = 0.01, to = 0.01),  # susceptance
        g = (from = 0.01, to = 0.01),  # conductance
        rating = 100.0,  # rating
        angle_limits = (min = -π/2, max = π/2),
    )

    line2 = PSY.Line(;
        name = "test_line_2",
        available = true,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = PSY.Arc(; from = bus1, to = bus2),
        r = 0.2,  # resistance
        x = 0.4,  # reactance
        b = (from = 0.02, to = 0.02),  # susceptance
        g = (from = 0.02, to = 0.02),  # conductance
        rating = 150.0,  # rating
        angle_limits = (min = -π/2, max = π/2),
    )

    # Create BranchesParallel
    bp = PNM.BranchesParallel([line1, line2])

    # Test get_equivalent_r and get_equivalent_x
    # For parallel branches: 1/Z_eq = 1/Z1 + 1/Z2
    # Z1 = 0.1 + 0.2j, Z2 = 0.2 + 0.4j
    # 1/Z1 = 1/(0.1+0.2j) = (0.1-0.2j)/((0.1)^2+(0.2)^2) = (0.1-0.2j)/0.05 = 2-4j
    # 1/Z2 = 1/(0.2+0.4j) = (0.2-0.4j)/((0.2)^2+(0.4)^2) = (0.2-0.4j)/0.2 = 1-2j
    # 1/Z_eq = (2-4j) + (1-2j) = 3-6j
    # Z_eq = 1/(3-6j) = (3+6j)/((3)^2+(6)^2) = (3+6j)/45 = 0.0667+0.1333j
    z1 = 0.1 + 0.2im
    z2 = 0.2 + 0.4im
    z_eq_expected = inv(inv(z1) + inv(z2))

    r_eq = PSY.get_r(bp)
    @test r_eq ≈ real(z_eq_expected) atol = 1e-6

    x_eq = PSY.get_x(bp)
    @test x_eq ≈ imag(z_eq_expected) atol = 1e-6

    # Test get_g: G_total = (from = G1_from + G2_from, to = G1_to + G2_to) = (from = 0.03, to = 0.03)
    g_eq = PSY.get_g(bp)
    @test g_eq.from ≈ 0.03 atol = 1e-6
    @test g_eq.to ≈ 0.03 atol = 1e-6

    # Test get_b: B_total = (from = B1_from + B2_from, to = B1_to + B2_to) = (from = 0.03, to = 0.03)
    b_eq = PSY.get_b(bp)
    @test b_eq.from ≈ 0.03 atol = 1e-6
    @test b_eq.to ≈ 0.03 atol = 1e-6

    # Test get_rating: Rating = (Rating1 + Rating2) / n = (100.0 + 150.0) / 2 = 125.0
    rating_eq = PSY.get_rating(bp)
    @test rating_eq ≈ 125.0 atol = 1e-6

    # Test get_available: all branches must be available
    @test PSY.get_available(bp) == true

    # Test get_α: average phase angle shift (both lines have 0.0)
    α_eq = PSY.get_α(bp)
    @test α_eq ≈ 0.0 atol = 1e-6

    # Clean up - don't add to system
    # (lines were created but not added to system)
end

@testset "Equivalent getters for BranchesSeries" begin
    # Create test system
    sys = PSB.build_system(PSB.PSITestSystems, "case10_radial_series_reductions")

    # Get buses for series configuration
    buses = collect(PSY.get_components(PSY.ACBus, sys))
    bus1 = buses[1]
    bus2 = buses[2]
    bus3 = buses[3]

    # Create test branches in series with specific values
    line1 = PSY.Line(;
        name = "series_line_1",
        available = true,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = PSY.Arc(; from = bus1, to = bus2),
        r = 0.05,  # resistance
        x = 0.1,   # reactance
        b = (from = 0.01, to = 0.0),  # susceptance
        g = (from = 0.01, to = 0.0),  # conductance
        rating = 100.0,  # rating
        angle_limits = (min = -π/2, max = π/2),
    )

    line2 = PSY.Line(;
        name = "series_line_2",
        available = true,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = PSY.Arc(; from = bus2, to = bus3),
        r = 0.15,  # resistance
        x = 0.3,   # reactance
        b = (from = 0.0, to = 0.02),  # susceptance
        g = (from = 0.0, to = 0.02),  # conductance
        rating = 80.0,  # rating
        angle_limits = (min = -π/2, max = π/2),
    )

    # Create BranchesSeries
    bs = PNM.BranchesSeries()
    PNM.add_branch!(bs, line1)
    PNM.add_branch!(bs, line2)

    # Test get_r: R_total = R1 + R2 = 0.05 + 0.15 = 0.2
    r_eq = PSY.get_r(bs)
    @test r_eq ≈ 0.2 atol = 1e-6

    # Test get_x: X_total = X1 + X2 = 0.1 + 0.3 = 0.4
    x_eq = PSY.get_x(bs)
    @test x_eq ≈ 0.4 atol = 1e-6

    # Test get_g: should be from endpoints (from = first.from, to = last.to)
    g_eq = PSY.get_g(bs)
    @test g_eq.from ≈ 0.01 atol = 1e-6
    @test g_eq.to ≈ 0.02 atol = 1e-6

    # Test get_b: should be from endpoints (from = first.from, to = last.to)
    b_eq = PSY.get_b(bs)
    @test b_eq.from ≈ 0.01 atol = 1e-6
    @test b_eq.to ≈ 0.02 atol = 1e-6

    # Test get_rating: Rating_total = min(Rating1, Rating2) = min(100.0, 80.0) = 80.0
    rating_eq = PSY.get_rating(bs)
    @test rating_eq ≈ 80.0 atol = 1e-6

    # Test get_available: all branches must be available
    @test PSY.get_available(bs) == true

    # Test get_α: sum of phase angle shifts (both lines have 0.0)
    α_eq = PSY.get_α(bs)
    @test α_eq ≈ 0.0 atol = 1e-6
end

@testset "Equivalent getters for ThreeWindingTransformerWinding" begin
    # Create a test system with three-winding transformers
    sys = PSB.build_system(PSB.PSITestSystems, "case10_radial_series_reductions")

    # Get a three-winding transformer from the system
    transformers_3w = collect(PSY.get_components(PSY.ThreeWindingTransformer, sys))

    if !isempty(transformers_3w)
        trf = first(transformers_3w)

        # Test winding 1 (primary)
        tw1 = PNM.ThreeWindingTransformerWinding(trf, 1)

        r1 = PSY.get_r(tw1)
        @test r1 == PSY.get_r_primary(trf)

        x1 = PSY.get_x(tw1)
        @test x1 == PSY.get_x_primary(trf)

        g1 = PSY.get_g(tw1)
        @test g1.from == real(PSY.get_g(trf))
        @test g1.to == 0.0

        b1 = PSY.get_b(tw1)
        @test b1.from == imag(PSY.get_g(trf))
        @test b1.to == 0.0

        rating1 = PSY.get_rating(tw1)
        # Should return winding-specific rating if non-zero, else transformer rating
        expected_rating1 = trf.rating_primary == 0.0 ? PSY.get_rating(trf) : trf.rating_primary
        @test rating1 == expected_rating1

        # Test get_available for winding 1
        @test PSY.get_available(tw1) == PSY.get_available(trf)

        # Test get_α for winding 1
        α1 = PSY.get_α(tw1)
        @test α1 == PSY.get_α_primary(trf)

        # Test winding 2 (secondary)
        tw2 = PNM.ThreeWindingTransformerWinding(trf, 2)

        r2 = PSY.get_r(tw2)
        @test r2 == PSY.get_r_secondary(trf)

        x2 = PSY.get_x(tw2)
        @test x2 == PSY.get_x_secondary(trf)

        g2 = PSY.get_g(tw2)
        @test g2.from == 0.0
        @test g2.to == 0.0

        b2 = PSY.get_b(tw2)
        @test b2.from == 0.0
        @test b2.to == 0.0

        rating2 = PSY.get_rating(tw2)
        # Should return winding-specific rating if non-zero, else transformer rating
        expected_rating2 = trf.rating_secondary == 0.0 ? PSY.get_rating(trf) : trf.rating_secondary
        @test rating2 == expected_rating2

        # Test get_α for winding 2
        α2 = PSY.get_α(tw2)
        @test α2 == PSY.get_α_secondary(trf)

        # Test winding 3 (tertiary)
        tw3 = PNM.ThreeWindingTransformerWinding(trf, 3)

        r3 = PSY.get_r(tw3)
        @test r3 == PSY.get_r_tertiary(trf)

        x3 = PSY.get_x(tw3)
        @test x3 == PSY.get_x_tertiary(trf)

        g3 = PSY.get_g(tw3)
        @test g3.from == 0.0
        @test g3.to == 0.0

        b3 = PSY.get_b(tw3)
        @test b3.from == 0.0
        @test b3.to == 0.0

        rating3 = PSY.get_rating(tw3)
        # Should return winding-specific rating if non-zero, else transformer rating
        expected_rating3 = trf.rating_tertiary == 0.0 ? PSY.get_rating(trf) : trf.rating_tertiary
        @test rating3 == expected_rating3

        # Test get_α for winding 3
        α3 = PSY.get_α(tw3)
        @test α3 == PSY.get_α_tertiary(trf)

        # Test error handling for invalid winding number
        tw_invalid = PNM.ThreeWindingTransformerWinding(trf, 4)
        @test_throws ArgumentError PSY.get_r(tw_invalid)
        @test_throws ArgumentError PSY.get_x(tw_invalid)
        @test_throws ArgumentError PSY.get_α(tw_invalid)
        # Note: get_rating and get_available don't throw for invalid winding as they use parent transformer values
    else
        @warn "No three-winding transformers found in test system; skipping some tests."
    end
end

@testset "Edge cases for equivalent getters" begin
    # Test with single branch in BranchesParallel
    sys = PSB.build_system(PSB.PSITestSystems, "case10_radial_series_reductions")
    bus1 = first(PSY.get_components(PSY.ACBus, sys))
    bus2 = collect(PSY.get_components(PSY.ACBus, sys))[2]

    line_single = PSY.Line(;
        name = "single_line",
        available = true,
        active_power_flow = 0.0,
        reactive_power_flow = 0.0,
        arc = PSY.Arc(; from = bus1, to = bus2),
        r = 0.1,
        x = 0.2,
        b = (from = 0.01, to = 0.01),
        g = (from = 0.01, to = 0.01),
        rating = 100.0,
        angle_limits = (min = -π/2, max = π/2),
    )

    bp_single = PNM.BranchesParallel([line_single])

    # For single branch, equivalent impedance calculation should still match branch values
    z_single = 0.1 + 0.2im
    z_eq_single = inv(inv(z_single))
    @test PSY.get_r(bp_single) ≈ real(z_eq_single) atol = 1e-6
    @test PSY.get_x(bp_single) ≈ imag(z_eq_single) atol = 1e-6
    @test PSY.get_g(bp_single).from ≈ 0.01 atol = 1e-6
    @test PSY.get_g(bp_single).to ≈ 0.01 atol = 1e-6
    @test PSY.get_b(bp_single).from ≈ 0.01 atol = 1e-6
    @test PSY.get_b(bp_single).to ≈ 0.01 atol = 1e-6
    @test PSY.get_rating(bp_single) ≈ 100.0 atol = 1e-6

    # Test with single branch in BranchesSeries
    bs_single = PNM.BranchesSeries()
    PNM.add_branch!(bs_single, line_single)

    @test PSY.get_r(bs_single) ≈ 0.1 atol = 1e-6
    @test PSY.get_x(bs_single) ≈ 0.2 atol = 1e-6
    @test PSY.get_g(bs_single).from ≈ 0.01 atol = 1e-6
    @test PSY.get_g(bs_single).to ≈ 0.01 atol = 1e-6
    @test PSY.get_b(bs_single).from ≈ 0.01 atol = 1e-6
    @test PSY.get_b(bs_single).to ≈ 0.01 atol = 1e-6
    @test PSY.get_rating(bs_single) ≈ 100.0 atol = 1e-6
end
