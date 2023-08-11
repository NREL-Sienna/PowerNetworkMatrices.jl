"""
Structure containing the BA matrix and other relevant data.

# Arguments
- `data::SparseArrays.SparseMatrixCSC{Float64, Int}`:
        the transposed BA matrix coming from the product between the Incidence 
        Matrix A and the Matrix of Susceptance B
- `axes<:NTuple{2, Dict}`:
        Tuple containing two vectors, the first one contains the names of each
        buse of the network (each one related to a row of the Matrix in "data"),
        the second one contains the names of each line of the network (each one
        related to a column of the Matrix in "data")
- `lookup<:NTuple{2, Dict}`:
        Tuple containing 2 Dictionaries mapping the number of rows and columns
        with the names of buses and branches
- `ref_bus_positions::Set{Int}`:
        Vector containing the indexes of the columns of the BA matrix corresponding
        to the refence buses
"""
struct BA_Matrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Float64}
    data::SparseArrays.SparseMatrixCSC{Float64, Int}
    axes::Ax
    lookup::L
    ref_bus_positions::Set{Int}
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
    axes = (bus_ax, line_ax)
    lookup = (make_ax_ref(bus_ax), make_ax_ref(line_ax))
    data = calculate_BA_matrix(branches, bus_lookup)
    return BA_Matrix(data, axes, lookup, ref_bus_positions)
end

"""
Structure containing the ABA matrix and other relevant data.

# Arguments
- `data::SparseArrays.SparseMatrixCSC{Float64, Int}`:
        the ABA matrix coming from the product between the Incidence Matrix A and
        the Matrix BA.
- `axes<:NTuple{2, Dict}`:
        Tuple containing two identical vectors, both containing the number of 
        each bus of the network (each one related to a row/column of the Matrix
        in "data"), excluding the slack buses.
- `lookup<:NTuple{2, Dict}`:
        Tuple containing 2 Dictionaries mapping the number of rows and columns
        with the number of the buses.
- `ref_bus_positions::Set{Int}`:
        Vector containing the indexes of the columns of the BA matrix corresponding
        to the refence buses
- `K<:Union{Nothing, KLU.KLUFactorization{Float64, Int}}`:
        either nothing or a container for KLU factorization matrices (LU factorization)
"""
struct ABA_Matrix{
    Ax,
    L <: NTuple{2, Dict},
    F <: Union{Nothing, KLU.KLUFactorization{Float64, Int}},
} <: PowerNetworkMatrix{Float64}
    data::SparseArrays.SparseMatrixCSC{Float64, Int}
    axes::Ax
    lookup::L
    ref_bus_positions::Set{Int}
    K::F
end

"""
Builds the ABA matrix from a System

# Arguments
- `sys::PSY.System`:
        system to consider

# Keyword arguments
- `factorize`: if true populates ABA_Matrix.K with KLU factorization matrices
"""
function ABA_Matrix(sys::PSY.System; factorize = false)
    branches = get_ac_branches(sys)
    buses = get_buses(sys)
    bus_lookup = make_ax_ref(buses)

    A, ref_bus_positions = calculate_A_matrix(branches, buses)
    BA = calculate_BA_matrix(branches, bus_lookup)
    ABA = calculate_ABA_matrix(A, BA, ref_bus_positions)

    bus_ax = [PSY.get_number(bus) for bus in buses]
    bus_ax_ = setdiff(bus_ax, ref_bus_positions)
    axes = (bus_ax_, bus_ax_)
    bus_ax_ref = make_ax_ref(bus_ax_)
    lookup = (bus_ax_ref, bus_ax_ref)
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

"""
Evaluates the LU factorization matrices of the ABA matrix, using KLU.

# Arguments
- `ABA::ABA_Matrix{Ax, L, Nothing} where {Ax, L <: NTuple{2, Dict}}`:
        container for the ABA matrix, with ABA.K == nothing (LU matrices in K not evaluated)
"""
function factorize(ABA::ABA_Matrix{Ax, L, Nothing}) where {Ax, L <: NTuple{2, Dict}}
    ABA_lu = ABA_Matrix(
        deepcopy(ABA.data),
        deepcopy(ABA.axes),
        deepcopy(ABA.lookup),
        deepcopy(ABA.ref_bus_positions),
        klu(ABA.data),
    )
    return ABA_lu
end

# checks if ABA has been factorized (if K contained LU matrices)
is_factorized(ABA::ABA_Matrix{Ax, L, Nothing}) where {Ax, L <: NTuple{2, Dict}} = false
is_factorized(
    ABA::ABA_Matrix{Ax, L, KLU.KLUFactorization{Float64, Int}},
) where {Ax, L <: NTuple{2, Dict}} = true

# get_index functions: BA_Matrix stores the transposed matrix, thus get index
# must export values according to [branch, bus] indexing.
function Base.getindex(A::BA_Matrix, line, bus)
    i, j = to_index(A, bus, line)
    return A.data[i, j]
end

function Base.getindex(
    A::BA_Matrix,
    line_number::Union{Int, Colon},
    bus_number::Union{Int, Colon},
)
    return A.data[bus_number, line_number]
end

# get_index functions: ABA_Matrix stores a square matrix whose number of rows 
# and column is equal to the number of the system's buses minus the slack ones,
# NOTE: bus_1, bus_2 are bus numbers not row and column indices!
function Base.getindex(A::ABA_Matrix, bus_1, bus_2)
    if bus_1 in A.ref_bus_positions || bus_2 in A.ref_bus_positions
        err_msg = string(" Rows and columns related to slack buses are not defined for the ABA matrix. \n",
            "Indices must be referred to any bus number that is not a slack one. \n",
            "For the current Systems the reference slack buses are: ", collect(A.ref_bus_positions), ". \n")
        error(err_msg)
    else
        i, j = to_index(A, bus_1, bus_2)
        return A.data[i, j]
    end
    return
end