module PowerNetworkMatrices

export ABA_Matrix
export BA_Matrix
export factorize
export find_subnetworks
export from_hdf5
export get_ptdf_data
export get_lodf_data
export get_bus_reduction_map
export get_radial_reduction
export get_ward_reduction
export get_reduction_type
export IncidenceMatrix
export is_factorized
export LODF
export PTDF
export NetworkReduction
export NetworkReductionTypes

export to_hdf5
export validate_connectivity
export VirtualLODF
export VirtualPTDF
export Ybus
export get_branch_lookups

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
include("network_reduction.jl")
include("ybus_calculations.jl")
include("incedence_matrix.jl")
include("breaker_switch_reduction.jl")
include("radial_reduction.jl")
include("degree_two_reduction.jl")
include("ward_reduction.jl")
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
