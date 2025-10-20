
@testset "LCCAdmittanceMatrix" begin
    sys = build_system(PSSEParsingTestSystems, "pti_two_terminal_hvdc_test_sys")
    lccs = get_components(TwoTerminalLCCLine, sys)
    lcc_admittance_matrix = PowerNetworkMatrices.LCCAdmittanceMatrix(collect(lccs), :FromTo)
    @test lcc_admittance_matrix[(1002, 1001), 1001] == ComplexF32(0.0)
    @test lcc_admittance_matrix[(1002, 1001), 1002] == ComplexF32(0.0)
    @test isa(lcc_admittance_matrix, PowerNetworkMatrices.LCCAdmittanceMatrix)
end
