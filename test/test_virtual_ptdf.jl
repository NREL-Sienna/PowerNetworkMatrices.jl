# @testset "Virtual PTDF matrices" begin
#     sys = PSB.build_system(PSB.PSYTestSystems, "tamu_ACTIVSg2000_sys")
#     ptdf_complete = PTDF(sys; linear_solver = "KLU")
#     ptdf_virtual = VirtualPTDF(sys)

#     for i in axes(ptdf_complete, 1)
#         comp = ptdf_complete[i, :]
#         virtual = ptdf_virtual[i, :]
#         for j in axes(ptdf_complete, 2)
#             # check values using PTDFs axes
#             @test isapprox(ptdf_complete[i, j], ptdf_virtual[i, j]; atol = 1e-10)
#         end
#     end
#     # Check the cache is populated
#     @test length(ptdf_virtual.cache) == length(ptdf_virtual.axes[1])
# end

## test new functions to consider different behavior #########################
using Revise
using Test
import Logging
using PowerNetworkMatrices
using PowerSystems
using InfrastructureSystems
using PowerSystemCaseBuilder
using TimeSeries
using BenchmarkTools
using ProfileView
using Cthulhu

const IS = InfrastructureSystems
const PSY = PowerSystems
const PSB = PowerSystemCaseBuilder

# sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")
sys5 = System("ACTIVSg2000.m")

# first test the original functions
ptdf_virtual = VirtualPTDF(sys5)

# check profiling now
ProfileView.@profview 
@code_warntype ptdf_virtual["TYLER 1 0-TYLER 2 0-i_3156", :]

empty!(ptdf_virtual.cache)
@code_warntype PowerNetworkMatrices._getindex(ptdf_virtual, 2, 2) 

"""
do test with solve also for the normal PTDF case to see if gains are present
"""