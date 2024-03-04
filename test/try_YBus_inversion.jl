using Revise
using Test
import Logging
import LinearAlgebra: I
using PowerNetworkMatrices
using PowerSystems
using InfrastructureSystems
using PowerSystemCaseBuilder
using TimeSeries
using KLU
using SparseArrays
using PowerFlows

const IS = InfrastructureSystems
const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder
const PNM = PowerNetworkMatrices

# ============================================================================
# ============================================================================

# try the code with a RTS system
sys = PSB.build_system(PSB.PSISystems, "RTS_GMLC_DA_sys");

# get network matrices
a_mat = PNM.IncidenceMatrix(sys)
ba_mat = PNM.BA_Matrix(sys)
ptdf_mat = PNM.PTDF(sys)

# get the 2T HVDC data, then the buses connecting
hvdc_line = PSY.get_component(PSY.TwoTerminalHVDCLine, sys, "DC1");
study_buses = [
    PSY.get_from(PSY.get_arc(hvdc_line)),
    PSY.get_to(PSY.get_arc(hvdc_line))
]

# get all other buses
other_buses = PSY.get_components(x -> !(x in study_buses), PSY.ACBus, sys)

# get flows and angles
data = PowerFlowData(DCPowerFlow(), sys)
power_injection =
    deepcopy(data.bus_activepower_injection - data.bus_activepower_withdrawals)
matrix_data = deepcopy(data.power_network_matrix.K)       # LU factorization of ABA
aux_network_matrix = deepcopy(data.aux_network_matrix)    # BA matrix

valid_ix = setdiff(1:length(power_injection), data.aux_network_matrix.ref_bus_positions)
ref_bus_angles = deepcopy(data.bus_angles)
ref_flow_values = deepcopy(data.branch_flow_values)

ref_bus_angles[valid_ix] = matrix_data \ power_injection[valid_ix]
ref_flow_values = transpose(aux_network_matrix.data) * ref_bus_angles

# get the angles from the PSY data
thetas = zeros((length(ba_mat.axes[1]),))
for b in PSY.get_components(PSY.ACBus, sys)
    idx = ba_mat.lookup[1][PSY.get_number(b)]
    thetas[idx] = PSY.get_angle(b)
end

# ! angles are different after computation


# now evaluate the Ward decomposition
ybus = PNM.Ybus(sys)
ward_data = Ward_data(ybus, PSY.get_number.(study_buses))
