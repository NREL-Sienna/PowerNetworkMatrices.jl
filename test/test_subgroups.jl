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

## check ABA matrix ##########################################################

# new version
ABA = ABA_Matrix(sys5)
# old version
BA = BA_Matrix(sys5)
A = IncidenceMatrix(sys5)
ABA_2 = PowerNetworkMatrices.calculate_ABA_matrix(A.data, BA.data, A.ref_bus_positions)
@test isapprox(ABA.data, ABA_2, atol=1e-7)

## compare M with ABA ########################################################

# get M matrix and test if the sub_networks are the same
M = AdjacencyMatrix(sys5)
sub1 = find_subnetworks(M)
sub2 = find_subnetworks(ABA)
# check that the methods have the same keys and values
@test keys(sub1) == keys(sub2)
for key in keys(sub1)
    @test setdiff(sub1[key], sub2[key]) == []
end

sys10 = get_10bus_test_system()

BA = BA_Matrix(sys5)

ABA = PowerNetworkMatrices.ABA_Matrix(sys5)

# get BA matrix
branches = PowerNetworkMatrices.get_ac_branches(sys5)
A = IncidenceMatrix(sys5)
BA = PowerNetworkMatrices.calculate_BA_matrix(
    branches, A.ref_bus_positions, A.lookup[2])

# that using first check 

# break down function
branches = PowerNetworkMatrices.get_ac_branches(sys5)
buses = PowerNetworkMatrices.get_buses(sys5)
ref_bus_positions = PowerNetworkMatrices.find_slack_positions(buses)
bus_lookup = PowerNetworkMatrices.make_ax_ref(buses)
line_ax = [PSY.get_name(branch) for branch in branches]
bus_ax = [PSY.get_number(bus) for bus in setdiff(buses,ref_bus_positions)]
axes = (line_ax, bus_ax)
data = PowerNetworkMatrices.calculate_BA_matrix(branches, ref_bus_positions, bus_lookup)


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