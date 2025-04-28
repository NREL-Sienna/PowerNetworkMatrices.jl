#= 
@testset "Breaker/Switch Reduction" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
    for line_name in ["1", "2"]
        l = get_component(Line, sys, line_name)
        set_r!(l, 0.0)
        set_x!(l, 0.0)
    end
    nr = get_breaker_switch_reduction(sys)
    @test get_bus_reduction_map(nr)[2] == Set([1, 4])
    @test get_removed_branches(nr) == Set(["1", "2"])
    @test get_reduction_type(nr) == [NetworkReductionTypes.BREAKER_SWITCH]
end

@testset "Separate reductions" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")
    l = get_component(ACBranch, sys, "Line12")
    set_r!(l, 0.0)
    set_x!(l, 0.0)
    nr1 = get_radial_reduction(sys)
    nr2 = get_breaker_switch_reduction(sys)
    nr3 = get_breaker_switch_reduction(sys; prior_reduction = nr1)
    @test nr1.retained_branches != nr2.retained_branches != nr3.retained_branches
    @test nr1.removed_branches != nr2.removed_branches != nr3.removed_branches
    isempty(nr3)
    nr4 = get_radial_reduction(sys; prior_reduction = nr2)
    @test nr3.retained_branches == nr4.retained_branches
    @test nr3.removed_branches == nr4.removed_branches
end

@testset "Overlapping reductions" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")
    l = get_component(ACBranch, sys, "Trans4")
    set_r!(l, 0.0)
    set_x!(l, 0.0)
    nr1 = get_radial_reduction(sys)
    nr2 = get_breaker_switch_reduction(sys)
    nr3 = get_breaker_switch_reduction(sys; prior_reduction = nr1)
    @test nr1.retained_branches == nr2.retained_branches == nr3.retained_branches
    @test nr1.removed_branches == nr2.removed_branches == nr3.removed_branches
end

@testset "Three reductions" begin
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")
    l = get_component(ACBranch, sys, "Line12")
    set_r!(l, 0.0)
    set_x!(l, 0.0)
    show_components(ACBus, sys, [:number])
    show_components(ACBus, sys, [:number])

    nr1 = get_breaker_switch_reduction(sys)
    nr2 = get_radial_reduction(sys; prior_reduction = nr1)
    nr3 = get_ward_reduction(sys, [7, 6, 12, 3, 1]; prior_reduction = nr2)  
    @test length(get_reduction_type(nr3)) == 3
end

@testset "Multi-reduction validation" begin 
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")
    nr1 = get_ward_reduction(sys,  [7, 6, 12, 3, 1])
    @test_throws IS.DataFormatError get_radial_reduction(sys; prior_reduction = nr1)
    @test_throws IS.DataFormatError get_breaker_switch_reduction(sys; prior_reduction = nr1)
    @test_throws IS.DataFormatError get_ward_reduction(sys, [1,2]; prior_reduction = nr1)
end  =#
