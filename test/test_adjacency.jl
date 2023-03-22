# TODO test M matrix with 10 bus case

# ! result: {1: [...]}, {2: [...]}

# TODO check consistency of dict values per key

# TODO use reference buses as key of the dictionary

using Revise
using Test
import Logging
using PowerNetworkMatrices
using PowerSystems
using InfrastructureSystems
using PowerSystemCaseBuilder
using TimeSeries

const IS = InfrastructureSystems
const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder

const BASE_DIR = dirname(dirname(Base.find_package("PowerSystems")))
const DATA_DIR = PSB.DATA_DIR

include("testing_data.jl")

# get systems
sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
sys10 = get_10bus_test_system()

# first check ptdf and M with the 5 bus case
M = AdjacencyMatrix(sys5)
ptdf = PTDF(sys5)
vptdf = VirtualPTDF(sys5)

# check subnetworks
ptdf.subnetworks
vptdf.subnetworks

# 10 bus case

ptdf10 = PTDF(sys10)
vptdf10 = VirtualPTDF(sys10)
a10 = IncidenceMatrix(sys10)
a10.data # * correct

PowerNetworkMatrices.get_ac_branches(sys10) # * correct

branches = PowerNetworkMatrices.get_ac_branches(sys10)
nodes = PowerNetworkMatrices.get_buses(sys10)
M, bus_ax_ref = PowerNetworkMatrices.calculate_adjacency(branches, nodes)

a10.ref_bus_positions
ptdf10.subnetworks
vptdf10.subnetworks

bus_numbers = [PSY.get_number(bus) for bus in nodes]
touched = []
bus_groups = Dict{Int, Vector{Int}}()

for i in axes(M, 1)
    if i ∉ touched
        # take bus numbers that are in the row
        cols = bus_numbers[M[i, :] .!= 0]
        if length(cols) == 1
            # islanded bus
            push!(touched, i)
            bus_groups[cols[1]] = cols
        elseif length(cols) > 1
            # subnetwork
            push!(touched, i)
            bus_groups[cols[1]] = cols
            for j in setdiff(axes(M, 1), touched)
                cols2 = bus_numbers[M[j, :] .!= 0]
                if length(cols2) > 1
                    if length(intersect(bus_groups[i], cols2)) != 0
                        push!(touched, j)
                        append!(bus_groups[i], setdiff(cols2, bus_groups[i]))
                    end
                end
            end
        end
    end
end
