# ! PROBLEMS:
# ! - vptsd[1, :] not working for calling or creating
# ! - vptsd["1", :] not working for creating

using Test
import Logging
using PowerNetworkMatrices
using PowerSystems
using InfrastructureSystems
using PowerSystemCaseBuilder
using TimeSeries
using ProfileView
using Cthulhu
using BenchmarkTools
using Profile

import PowerNetworkMatrices: get_ac_branches
import PowerNetworkMatrices: get_buses, make_ax_ref

const IS = InfrastructureSystems
const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder

# break down the code

## VIRTUAL PTDF CREATION #####################################################
sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")

# get branches and buses
branches = get_ac_branches(sys5)
nodes = get_buses(sys5)

# get virtual PTDF

dist_slack::Vector{Float64} = [0.1]
tol::Float64 = eps()

# new
@btime begin
    line_ax = String[PSY.get_name(branch) for branch in branches]
    bus_ax = Int64[PSY.get_number(bus) for bus in nodes]
    axes = Tuple{Vector{String}, Vector{Int64}}[(line_ax, bus_ax)]
    look_up = Tuple{Dict{String, Int64}, Dict{Int64, Int64}}[
        (make_ax_ref(line_ax), make_ax_ref(bus_ax))]
end
"""
4.101 μs (27 allocations: 1.83 KiB)
"""

ProfileView.@profview begin
    line_ax = String[PSY.get_name(branch) for branch in branches]
    bus_ax = Int64[PSY.get_number(bus) for bus in nodes]
    axes = Tuple{Vector{String}, Vector{Int64}}[(line_ax, bus_ax)]
    look_up = Tuple{Dict{String, Int64}, Dict{Int64, Int64}}[
        (make_ax_ref(line_ax), make_ax_ref(bus_ax))]
end

# old
@btime begin
    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in nodes]
    axes = (line_ax, bus_ax)
    look_up = (make_ax_ref(line_ax), make_ax_ref(bus_ax))
end
"""
1.188 μs (15 allocations: 1.33 KiB)
old is better despite worse profiling plot
"""

ProfileView.@profview begin
    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in nodes]
    axes = (line_ax, bus_ax)
    look_up = (make_ax_ref(line_ax), make_ax_ref(bus_ax))
end

ProfileView.@profview begin
    A, slack_positions = PowerNetworkMatrices.calculate_A_matrix(branches, nodes)
    BA = PowerNetworkMatrices.calculate_BA_matrix(branches, slack_positions, make_ax_ref(nodes))
    ABA = PowerNetworkMatrices.calculate_ABA_matrix(A, BA, slack_positions)
    empty_cache = Dict{Int, Array{Float64}}()
    temp_data = Vector{Float64}(undef, size(BA, 2))
end
"""
this is ok, however can't understand why most of the time is spent out of the functions
"""

# try now getindex method
ptdf_virtual = VirtualPTDF(branches, nodes)
# evaluate new value
ProfileView.@profview ptdf_virtual[1, 1]
# get stored value
Profile.init()
ProfileView.@profview ptdf_virtual[1, 1]
