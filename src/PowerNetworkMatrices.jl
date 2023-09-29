module PowerNetworkMatrices

export ABA_Matrix
export AdjacencyMatrix
export BA_Matrix
export factorize
export find_subnetworks
export from_hdf5
export get_ptdf_data
export get_lodf_data
export IncidenceMatrix
export is_factorized
export LODF
export PTDF
export to_hdf5
export validate_connectivity
export VirtualLODF
export VirtualPTDF
export Ybus

using DocStringExtensions
#import MKL
import InfrastructureSystems
import PowerSystems
import PowerSystems: ACBusTypes

const IS = InfrastructureSystems
const PSY = PowerSystems

import SparseArrays
import SparseArrays: rowvals, nzrange
import HDF5
import KLU: klu
import KLU
import LinearAlgebra
import LinearAlgebra: BLAS.gemm
import LinearAlgebra: ldiv!, mul!, I, dot
import Pardiso

@template DEFAULT = """
                    $(SIGNATURES)
                    $(DOCSTRING)
                    """

# network calculations
include("PowerNetworkMatrix.jl")
include("BA_ABA_matrices.jl")
include("incedence_matrix.jl")
include("adjacency_matrix.jl")
include("common.jl")
include("definitions.jl")
include("ptdf_calculations.jl")
include("ybus_calculations.jl")
include("row_cache.jl")
include("virtual_ptdf_calculations.jl")
include("lodf_calculations.jl")
include("virtual_lodf_calculations.jl")
include("system_utils.jl")
include("serialization.jl")

end
