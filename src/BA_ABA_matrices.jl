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
    subnetwork_axes::Dict{Int, Ax}
    network_reduction_data::NetworkReductionData
end

get_axes(M::BA_Matrix) = M.axes
get_lookup(M::BA_Matrix) = M.lookup
get_ref_bus(M::BA_Matrix) = sort!(collect(keys(M.subnetwork_axes)))
get_ref_bus_position(M::BA_Matrix) = [get_bus_lookup(M)[x] for x in keys(M.subnetwork_axes)]
get_network_reduction_data(M::BA_Matrix) = M.network_reduction_data
get_bus_axis(M::BA_Matrix) = M.axes[1]
get_bus_lookup(M::BA_Matrix) = M.lookup[1]
get_arc_axis(M::BA_Matrix) = M.axes[2]
get_arc_lookup(M::BA_Matrix) = M.lookup[2]
stores_transpose(::BA_Matrix) = true

function BA_Matrix(sys::PSY.System;
    network_reductions::Vector{NetworkReduction} = Vector{NetworkReduction}(),
    kwargs...,
)
    return BA_Matrix(
        Ybus(
            sys;
            network_reductions = network_reductions,
            kwargs...,
        ),
    )
end

function BA_Matrix(ybus::Ybus)
    nr = ybus.network_reduction_data
    bus_ax = get_bus_axis(ybus)
    bus_lookup = get_bus_lookup(ybus)
    arc_ax = get_arc_axis(nr)
    n_isolated_buses = length(get_isolated_buses(ybus))
    n_entries = length(arc_ax) * 2 + n_isolated_buses
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
    BA_I[(end - n_isolated_buses + 1):end] =
        [get_bus_index(x, bus_lookup, nr) for x in get_isolated_buses(ybus)]
    BA_J[(end - n_isolated_buses + 1):end] = ones(n_isolated_buses)
    BA_V[(end - n_isolated_buses + 1):end] = zeros(n_isolated_buses)
    data = SparseArrays.sparse(BA_I, BA_J, BA_V)
    SparseArrays.dropzeros!(data)
    axes = (bus_ax, arc_ax)
    lookup = (make_ax_ref(bus_ax), make_ax_ref(arc_ax))
    subnetwork_axes = make_bus_arc_subnetwork_axes(ybus)
    return BA_Matrix(data, axes, lookup, subnetwork_axes, ybus.network_reduction_data)
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
    subnetwork_axes::Dict{Int, Ax}
    ref_bus_position::Vector{Int}
    K::F
    network_reduction_data::NetworkReductionData
end

get_axes(M::ABA_Matrix) = M.axes
get_lookup(M::ABA_Matrix) = M.lookup
get_ref_bus(M::ABA_Matrix) = sort!(collect(keys(M.subnetwork_axes)))
get_ref_bus_position(M::ABA_Matrix) = M.ref_bus_position
get_network_reduction_data(M::ABA_Matrix) = M.network_reduction_data
get_bus_axis(M::ABA_Matrix) = M.axes[1]
get_bus_lookup(M::ABA_Matrix) = M.lookup[1]

function ABA_Matrix(sys::PSY.System;
    factorize::Bool = false,
    network_reductions::Vector{NetworkReduction} = NetworkReduction[],
    kwargs...,
)
    ymatrix = Ybus(
        sys;
        network_reductions = network_reductions,
        kwargs...,
    )
    ref_bus_positions = get_ref_bus_position(ymatrix)
    A = IncidenceMatrix(ymatrix)
    BA = BA_Matrix(ymatrix)
    ABA = calculate_ABA_matrix(A.data, BA.data, Set(ref_bus_positions))
    axes, subnetwork_axes = _remake_axes_without_ref(ymatrix.axes, ymatrix.subnetwork_axes)
    bus_ax_ref = make_ax_ref(axes[1])
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
        subnetwork_axes,
        ref_bus_positions,
        K,
        ymatrix.network_reduction_data,
    )
end

function _remake_axes_without_ref(
    axes::Tuple{Vector{Int}, Vector{Int}},
    subnetwork_axes::Dict{Int, Tuple{Vector{Int}, Vector{Int}}},
)
    ref_bus_numbers = collect(keys(subnetwork_axes))
    bus_ax_ = setdiff(axes[1], Set(ref_bus_numbers))
    axes_ref = (bus_ax_, bus_ax_)
    subnetwork_axes_ref = deepcopy(subnetwork_axes)
    for (k, v) in subnetwork_axes_ref
        setdiff!(v[1], [k])
        setdiff!(v[2], [k])
    end
    return axes_ref, subnetwork_axes_ref
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
        deepcopy(ABA.subnetwork_axes),
        deepcopy(ABA.ref_bus_position),
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
    ref_bus_numbers = Set(get_ref_bus(A))
    if _is_ref_bus(bus_1, ref_bus_numbers) || _is_ref_bus(bus_2, ref_bus_numbers)
        err_msg = string(
            " Rows and columns related to slack buses are not defined for the ABA matrix. \n",
            "Indices must be referred to any bus number that is not a slack one. \n",
            "For the current Systems the reference slack buses are: ",
            collect(ref_bus_numbers), ". \n")
        error(err_msg)
    else
        i, j = to_index(A, bus_1, bus_2)
        return A.data[i, j]
    end
    return
end
