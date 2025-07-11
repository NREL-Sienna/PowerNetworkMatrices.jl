"""
Structure containing the BA matrix and other relevant data.

# Arguments
- `data::SparseArrays.SparseMatrixCSC{Float64, Int}`:
        the transposed BA matrix coming from the product between the Incidence
        Matrix A and the Matrix of Susceptance B
- `axes<:NTuple{2, Dict}`:
        Tuple containing two vectors, the first one contains the names of each
        buses of the network (each one related to a row of the Matrix in "data"),
        the second one contains the names of each line of the network (each one
        related to a column of the Matrix in "data")
- `lookup<:NTuple{2, Dict}`:
        Tuple containing 2 Dictionaries mapping the number of rows and columns
        with the names of buses and branches
- `ref_bus_positions::Set{Int}`:
        Set containing the indexes of the columns of the BA matrix corresponding
        to the reference buses
- `network_reduction::NetworkReduction`:
        Structure containing the details of the network reduction applied when computing the matrix
"""
struct BA_Matrix{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Float64}
    data::SparseArrays.SparseMatrixCSC{Float64, Int}
    axes::Ax
    lookup::L
    ref_bus_positions::Set{Int}
    network_reduction_data::NetworkReductionData
end

function BA_Matrix(sys::PSY.System;
    check_connectivity::Bool = true,
    network_reductions::Vector{NetworkReduction} = Vector{NetworkReduction}(),
    kwargs...,
)
    return BA_Matrix(
        Ybus(
            sys;
            check_connectivity = check_connectivity,
            network_reductions = network_reductions,
            kwargs...,
        ),
    )
end

function BA_Matrix(ybus::Ybus)
    nr = ybus.network_reduction_data
    direct_arcs = [x for x in keys(nr.direct_branch_map)]
    parallel_arcs = [x for x in keys(nr.parallel_branch_map)]
    series_arcs = [x for x in keys(nr.series_branch_map)]
    transformer_arcs = [x for x in keys(nr.transformer3W_map)]
    bus_ax = ybus.axes[1]
    bus_lookup = ybus.lookup[1]
    arc_ax = vcat(direct_arcs, parallel_arcs, series_arcs, transformer_arcs)
    n_entries = length(arc_ax) * 2
    BA_I = Vector{Int}(undef, n_entries)
    BA_J = Vector{Int}(undef, n_entries)
    BA_V = Vector{Float64}(undef, n_entries)
    for (ix_arc, arc) in enumerate(arc_ax)
        ix_from_bus = get_bus_index(arc[1], bus_lookup, nr)
        ix_to_bus = get_bus_index(arc[2], bus_lookup, nr)
        Yt = -1 * ybus.data[ix_from_bus, ix_to_bus]
        Zt = 1 / Yt
        # TODO - should we consider phase shift?
        b = 1 / imag(Zt)
        BA_I[2 * ix_arc - 1] = ix_from_bus
        BA_J[2 * ix_arc - 1] = ix_arc
        BA_V[2 * ix_arc - 1] = b
        BA_I[2 * ix_arc] = ix_to_bus
        BA_J[2 * ix_arc] = ix_arc
        BA_V[2 * ix_arc] = -1 * b
    end
    data = SparseArrays.sparse(BA_I, BA_J, BA_V)
    axes = (bus_ax, arc_ax)
    lookup = (make_ax_ref(bus_ax), make_ax_ref(arc_ax))
    ref_bus_positions = Set([lookup[1][x] for x in ybus.ref_bus_numbers])
    return BA_Matrix(data, axes, lookup, ref_bus_positions, ybus.network_reduction_data)
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
        to the reference buses
- `K<:Union{Nothing, KLU.KLUFactorization{Float64, Int}}`:
        either nothing or a container for KLU factorization matrices (LU factorization)
- `network_reduction::NetworkReduction`:
        Structure containing the details of the network reduction applied when computing the matrix
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
    ref_bus_numbers::Set{Int}
    K::F
    network_reduction_data::NetworkReductionData
end

function ABA_Matrix(sys::PSY.System;
    factorize::Bool = false,
    check_connectivity::Bool = true,
    network_reductions::Vector{NetworkReduction} = NetworkReduction[],
    kwargs...,
)
    ymatrix = Ybus(
        sys;
        check_connectivity = check_connectivity,
        network_reductions = network_reductions,
        kwargs...,
    )
    ref_bus_positions = Set([ymatrix.lookup[1][x] for x in ymatrix.ref_bus_numbers])
    A = IncidenceMatrix(ymatrix)
    BA = BA_Matrix(ymatrix)
    ABA = calculate_ABA_matrix(A.data, BA.data, ref_bus_positions)
    bus_ax = ymatrix.axes[1]
    ref_bus_numbers = ymatrix.ref_bus_numbers
    bus_ax_ = setdiff(bus_ax, ref_bus_numbers)
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
        ref_bus_numbers,
        K,
        ymatrix.network_reduction_data,
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
        deepcopy(ABA.ref_bus_numbers),
        klu(ABA.data),
        deepcopy(ABA.network_reduction_data),
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

_is_ref_bus(bus::Int, ref_bus_numbers::Set{Int}) = (bus in ref_bus_numbers)
_is_ref_bus(bus::PSY.ACBus, ref_bus_numbers::Set{Int}) =
    (PSY.get_number(bus) in ref_bus_numbers)
_is_ref_bus(::Colon, ::Set{Int}) = false
# get_index functions: ABA_Matrix stores a square matrix whose number of rows
# and column is equal to the number of the system's buses minus the slack ones,
# NOTE: bus_1, bus_2 are bus numbers (or PSY.ACBus objects) not row and column indices!
function Base.getindex(A::ABA_Matrix, bus_1, bus_2)
    if _is_ref_bus(bus_1, A.ref_bus_numbers) || _is_ref_bus(bus_2, A.ref_bus_numbers)
        err_msg = string(
            " Rows and columns related to slack buses are not defined for the ABA matrix. \n",
            "Indices must be referred to any bus number that is not a slack one. \n",
            "For the current Systems the reference slack buses are: ",
            collect(A.ref_bus_numbers), ". \n")
        error(err_msg)
    else
        i, j = to_index(A, bus_1, bus_2)
        return A.data[i, j]
    end
    return
end
