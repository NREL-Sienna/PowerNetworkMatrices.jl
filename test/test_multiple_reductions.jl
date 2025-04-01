
@testset "Zero Impedance Reduction" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    for line_name in ["1", "2"]
        l = get_component(Line, sys, line_name)
        set_r!(l, 0.0)
        set_x!(l, 0.0)
    end
    zr = get_zero_impedance_reduction(sys)
    @test get_bus_reduction_map(zr)[2] == Set([1, 4])
    @test get_removed_branches(zr) == Set(["1", "2"])
    @test get_reduction_type(zr) == [NetworkReductionTypes.ZERO_IMPEDANCE]
end
