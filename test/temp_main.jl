using PowerSystemCaseBuilder
using PowerSystems

const PNM = include("/Users/acastel2/Documents/GitHub/PowerNetworkMatrices.jl/src/PowerNetworkMatrices.jl")
const PSY = PowerSystems


# get system
sys = build_system(PSITestSystems, "test_RTS_GMLC_sys")

branches = PNM.get_ac_branches(sys)
nodes = PNM.get_buses(sys)
S, A = PNM._buildptdf(branches, nodes)
ptdf, a = PNM.PTDF(sys)

A = PNM.IncidenceMatrix(sys)

BA = PNM.calculate_BA_matrix(PNM.get_ac_branches(sys),
                             PNM.get_slack_position(A),
                             PNM.get_lookup(A))
ABA = PNM.calculate_ABA_matrix(A.data, BA)


branches = get_ac_branches(sys)
buses = get_buses(sys)

PSY.get_number.(buses)

@code_warntype PNM.get_ac_branches(sys)
@code_warntype PNM.get_buses(sys)

