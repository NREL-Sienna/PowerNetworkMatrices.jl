@testset "Test 3W Integration to Ybus Matrix" begin
    sys = build_system(PSSEParsingTestSystems, "pti_three_winding_test_2_sys")

    bus2 = get_component(ACBus, sys, "FAV SPOT 02")

    fixed_y = FixedAdmittance(;
        name = "fixed",
        available = true,
        bus = bus2,
        Y = 0.0 + 10im,
    )
    add_component!(sys, fixed_y)

    switch_y = SwitchedAdmittance(;
        name = "switch",
        available = true,
        bus = bus2,
        Y = 0.0 + 4im,
        number_of_steps = [2],
        Y_increase = [0.0 + 0.1im],
    )
    add_component!(sys, switch_y)

    ybus = Ybus(sys)

    @test isapprox(ybus[1001, 1001], 9.0e-9 - 20.0im, atol = 1e-4)
end
