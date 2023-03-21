# TODO test M matrix with 10 bus case

# ! result: {1: [...]}, {2: [...]}

# TODO check consistency of dict values per key

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
