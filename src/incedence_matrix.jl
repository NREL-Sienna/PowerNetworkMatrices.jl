"""
Power Transfer Distribution Factors (PTDF) indicate the incremental change in real power that occurs on transmission lines due to real power injections changes at the buses.

The A struct is indexed using the Bus numbers and branch names
"""
struct IncidenceMatrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{T <: Int}
    data::Array{Matrix, 2}
    axes::Ax
    lookup::L
end

function IncidenceMatrix(sys::PSY.System)
    call methods in common.jl

    data = calculate_A_matrix(branches, buses)
    return IncidenceMatrix(data, axes, lookup)
end
