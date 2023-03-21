using Revise
using Test
import Logging
using PowerNetworkMatrices
using PowerSystems
using InfrastructureSystems
using PowerSystemCaseBuilder
using TimeSeries
using SuiteSparse

const IS = InfrastructureSystems
const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder

const BASE_DIR = dirname(dirname(Base.find_package("PowerSystems")))
const DATA_DIR = PSB.DATA_DIR

include("testing_data.jl")

sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
sys10 = get_10bus_test_system()

M = AdjacencyMatrix(sys5)
BA = PowerNetworkMatrices.BA_matrix(sys5)

# break down function
branches = PowerNetworkMatrices.get_ac_branches(sys5)
buses = PowerNetworkMatrices.get_buses(sys5)
ref_bus_positions = PowerNetworkMatrices.find_slack_positions(buses)
bus_lookup = PowerNetworkMatrices.make_ax_ref(buses)
line_ax = [PSY.get_name(branch) for branch in branches]
bus_ax = [PSY.get_number(bus) for bus in setdiff(buses,ref_bus_positions)]
axes = (line_ax, bus_ax)
data = PowerNetworkMatrices.calculate_BA_matrix(branches, ref_bus_positions, bus_lookup)

typeof(axes)

rows = SparseArrays.rowvals(M)
_, n = size(M)
touched = Set{Int}()
bus_groups = Dict{Int, Set{Int}}()

for j in 1:n
    row_ix = rows[j]
    if bus_numbers[row_ix] ∉ touched
        push!(touched, bus_numbers[row_ix])
        bus_groups[bus_numbers[row_ix]] = Set{Int}()
        dfs(row_ix, M, bus_numbers, bus_groups[bus_numbers[row_ix]], touched)
    end
end


# change key to reference bus

for key in keys(bus_groups)
    
end

bus_groups