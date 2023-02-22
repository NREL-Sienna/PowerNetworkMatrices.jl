module PowerNetworkMatrices

export PTDF
export Ybus
export LODF
export Adjacency

export dfs_connectivity
export find_connected_components
export validate_connectivity

using DocStringExtensions
import InfrastructureSystems
import PowerSystems
import PowerSystems: BusTypes

const IS = InfrastructureSystems
const PSY = PowerSystems

import SparseArrays
import LinearAlgebra: LAPACK.getri!
import LinearAlgebra: LAPACK.getrf!
import LinearAlgebra: BLAS.gemm
import LinearAlgebra

@template (FUNCTIONS, METHODS) = """
                                 $(TYPEDSIGNATURES)
                                 $(DOCSTRING)
                                 """

# network calculations
include("PowerNetworkMatrix.jl")
include("common.jl")
include("ybus_calculations.jl")
include("ptdf_calculations.jl")
include("lodf_calculations.jl")

end
