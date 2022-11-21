module PowerNetworkMatrices

export PTDF
export Ybus
export LODF
export Adjacency

export validate_connectivity
export find_connected_components

using DocStringExtensions
import PowerSystems
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
include("common.jl")
include("ybus_calculations.jl")
include("ptdf_calculations.jl")
include("lodf_calculations.jl")

end
