"""
Branch admittance matrix

# Arguments
- `data::SparseArrays.SparseMatrixCSC{ComplexF32, Int}`:
        The branch admittance matrix in the from-to direction
- `axes<:NTuple{2, Dict}`:
        Tuple containing two vectors (the first one showing the arc tuples,
        the second showing the buses numbers).
- `lookup<:NTuple{2, Dict}`:
        Tuple containing two dictionaries, the first mapping the arc tuples
        and the second the buses with their enumerated indexes.
- `network_reduction::NetworkReduction`:
        Structure containing the details of the network reduction applied when computing the matrix
- `direction::Symbol`:
        Direction of admittance (:FromTo or :ToFrom)
"""
struct BranchAdmittanceMatrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{ComplexF32}
    data::SparseArrays.SparseMatrixCSC{ComplexF32, Int}
    axes::Ax
    lookup::L
    network_reduction_data::NetworkReductionData
    direction::Symbol
end

# functions to get stored data
get_axes(M::BranchAdmittanceMatrix) = M.axes
get_lookup(M::BranchAdmittanceMatrix) = M.lookup
get_ref_bus(M::BranchAdmittanceMatrix) = sort!(collect(keys(M.subnetwork_axes)))
get_ref_bus_position(M::BranchAdmittanceMatrix) =
    [get_bus_lookup(M)[x] for x in keys(M.subnetwork_axes)]
get_network_reduction_data(M::BranchAdmittanceMatrix) = M.network_reduction_data
get_arc_axis(M::BranchAdmittanceMatrix) = M.axes[1]
get_arc_lookup(M::BranchAdmittanceMatrix) = M.lookup[1]
get_bus_axis(M::BranchAdmittanceMatrix) = M.axes[2]
get_bus_lookup(M::BranchAdmittanceMatrix) = M.lookup[2]
get_direction(M::BranchAdmittanceMatrix) = M.direction

# TODO - define interface for getting branch admittance matrices from system
