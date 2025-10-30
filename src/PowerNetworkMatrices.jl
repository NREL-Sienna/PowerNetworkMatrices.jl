module PowerNetworkMatrices

export ABA_Matrix
export BA_Matrix
export factorize
export find_subnetworks
export from_hdf5
export get_ptdf_data
export get_lodf_data
export get_bus_reduction_map
export get_ward_reduction
export get_reductions
export get_network_reduction_data
export IncidenceMatrix
export AdjacencyMatrix
export is_factorized
export LODF
export PTDF
export NetworkReduction
export NetworkReductionData
export RadialReduction
export DegreeTwoReduction
export WardReduction

export depth_first_search
export iterative_union_find
export to_hdf5
export validate_connectivity
export VirtualLODF
export VirtualPTDF
export Ybus
export ArcAdmittanceMatrix

using DocStringExtensions
import InfrastructureSystems
import PowerSystems
import PowerSystems: ACBusTypes

const IS = InfrastructureSystems
const PSY = PowerSystems

@static if (Sys.ARCH === :x86_64 || Sys.ARCH === :i686) && !Sys.isapple()
    using MKL
    using Pardiso
    const USE_MKL = MKL.MKL_jll.is_available()
else
    const USE_MKL = false
end

@static if Sys.isapple()
    using AppleAccelerate
    const USE_AA = true
else
    const USE_AA = false
end

import SparseArrays
import SparseArrays: rowvals, nzrange
import HDF5
import KLU: klu
import KLU
import LinearAlgebra
import LinearAlgebra: BLAS.gemm
import LinearAlgebra: ldiv!, mul!, I, dot
import LinearAlgebra: LAPACK.getrf!, LAPACK.getrs!

@template DEFAULT = """
                    $(SIGNATURES)
                    $(DOCSTRING)
                    """

# network calculations
include("PowerNetworkMatrix.jl")
include("definitions.jl")
include("ThreeWindingTransformerWinding.jl")
include("NetworkReduction.jl")
include("radial_reduction.jl")
include("degree_two_reduction.jl")
include("ward_reduction.jl")
include("ReductionContainer.jl")
include("NetworkReductionData.jl")
include("ArcAdmittanceMatrix.jl")
include("Ybus.jl")
include("IncidenceMatrix.jl")
include("AdjacencyMatrix.jl")
include("connectivity_checks.jl")
include("subnetworks.jl")
include("common.jl")
include("BA_ABA_matrices.jl")
include("ptdf_calculations.jl")
include("row_cache.jl")
include("virtual_ptdf_calculations.jl")
include("lodf_calculations.jl")
include("virtual_lodf_calculations.jl")
include("system_utils.jl")
include("serialization.jl")

end
