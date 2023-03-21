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

## 5 bus case ################################################################

sys5 = PSB.build_system(PSB.PSITestSystems, "c_sys5")


## compare M with ABA ########################################################

# get M matrix and test if the sub_networks are the same
M = AdjacencyMatrix(sys5)
sub1 = find_subnetworks(M)

# TODO: add test for 2 network
# TODO: keys of the dictionary must be the reference buses

## 10 bus case ################################################################

sys10 = get_10bus_test_system()