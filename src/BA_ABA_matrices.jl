"""
Structure containing the BA matrix and other relevant data.

# Arguments
- `data::SparseArrays.SparseMatrixCSC{Float64, Int}`:
        the BA matrix coming from the product between the Incidence Matrix A and
        the Matrix of Susceptance B
- `axes<:NTuple{2, Dict}`:
        Tuple containing two vectors, the first one contains the names of each 
        line of the network (each one related to a row of the Matrix in "data"),
        the second one contains the names of each bus of the network (each one
        related to a column of the Matrix in "data")
- `lookup<:NTuple{2, Dict}`:
        Tuple containing 2 Dictionaries mapping the number of rows and columns 
        with the names of branches and buses
- `ref_bus_positions::Vector{Int}`:
        Vector containing the indexes of the columns of the BA matrix corresponding
        to the refence buses
"""
struct BA_Matrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Float64}
    data::SparseArrays.SparseMatrixCSC{Float64, Int}
    axes::Ax
    lookup::L
    ref_bus_positions::Vector{Int}
end

"""
Build the BA matrix from a given System
"""
function BA_Matrix(sys::PSY.System)
    branches = get_ac_branches(sys)
    buses = get_buses(sys)
    ref_bus_positions = find_slack_positions(buses)
    bus_lookup = make_ax_ref(buses)
    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in setdiff(buses, ref_bus_positions)]
    axes = (line_ax, bus_ax)
    lookup = (make_ax_ref(line_ax), make_ax_ref(bus_ax))
    data = calculate_BA_matrix(branches, ref_bus_positions, bus_lookup)
    return BA_Matrix(data, axes, lookup, ref_bus_positions)
end

"""
Structure containing the ABA matrix and other relevant data.

# Arguments
- `data::SparseArrays.SparseMatrixCSC{Float64, Int}`:
        the ABA matrix coming from the product between the Incidence Matrix A and
        the Matrix BA
- `axes<:NTuple{2, Dict}`:
        Tuple containing two vectors, the first one contains the names of each 
        line of the network (each one related to a row of the Matrix in "data"),
        the second one contains the names of each bus of the network (each one
        related to a column of the Matrix in "data")
- `lookup<:NTuple{2, Dict}`:
        Tuple containing 2 Dictionaries mapping the number of rows and columns 
        with the names of branches and buses
- `ref_bus_positions::Vector{Int}`:
        Vector containing the indexes of the columns of the BA matrix corresponding
        to the refence buses
"""
struct ABA_Matrix{
    Ax,
    L <: NTuple{2, Dict},
    F <: Union{Nothing, KLU.KLUFactorization{Float64, Int}},
} <: PowerNetworkMatrix{Float64}
    data::SparseArrays.SparseMatrixCSC{Float64, Int}
    axes::Ax
    lookup::L
    ref_bus_positions::Vector{Int}
    K::F
end

"""
Builds the ABA matrix from a System
"""
function ABA_Matrix(sys::PSY.System; factorize = false)
    branches = get_ac_branches(sys)
    buses = get_buses(sys)
    bus_lookup = make_ax_ref(buses)

    A, ref_bus_positions = calculate_A_matrix(branches, buses)
    BA = calculate_BA_matrix(branches, ref_bus_positions, bus_lookup)
    ABA = calculate_ABA_matrix(A, BA, ref_bus_positions)

    line_ax = [PSY.get_name(branch) for branch in branches]
    bus_ax = [PSY.get_number(bus) for bus in buses]
    axes = (line_ax, setdiff(bus_ax, ref_bus_positions))
    lookup = (make_ax_ref(line_ax), make_ax_ref(bus_ax))
    if factorize
        K = klu(ABA)
    else
        K = nothing
    end
    return ABA_Matrix(
        ABA,
        axes,
        lookup,
        ref_bus_positions,
        K,
    )
end

function factorize(ABA::ABA_Matrix{Ax, L, Nothing}) where {Ax, L <: NTuple{2, Dict}}
    return ABA_Matrix(
        deepcopy(ABA.data),
        deepcopy(ABA.axes),
        deepcopy(ABA.lookup),
        deepcopy(ABA.ref_bus_positions),
        KLU(ABA.data),
    )
end

is_factorized(ABA::ABA_Matrix{Ax, L, Nothing}) where {Ax, L <: NTuple{2, Dict}} = false
is_factorized(
    ABA::ABA_Matrix{Ax, L, KLU.KLUFactorization{Float64, Int}},
) where {Ax, L <: NTuple{2, Dict}} = true
