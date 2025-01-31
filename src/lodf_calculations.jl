"""
The Line Outage Distribution Factor (LODF) matrix gathers a sensitivity coefficients
of how a change in a line's flow affects the flows on other lines in the system.

# Arguments
- `data<:AbstractArray{Float64, 2}`:
        the transposed LODF matrix.
- `axes<:NTuple{2, Dict}`:
        Tuple containing two identical vectors containing the names of the
        branches related to each row/column.
- `lookup<:NTuple{2, Dict}`:
        Tuple containing two identical dictionaries mapping the branches
         their enumerated indexes (row and column numbers).
- `tol::Base.RefValue{Float64}`:
        tolerance used for sparsifying the matrix (dropping element whose
        absolute value is below this threshold).
- `network_reduction::NetworkReduction`:
        Structure containing the details of the network reduction applied when computing the matrix
"""
struct LODF{Ax, L <: NTuple{2, Dict}, M <: AbstractArray{Float64, 2}} <:
       PowerNetworkMatrix{Float64}
    data::M
    axes::Ax
    lookup::L
    tol::Base.RefValue{Float64}
    network_reduction::NetworkReduction
end

function _buildlodf(
    a::SparseArrays.SparseMatrixCSC{Int8, Int},
    ptdf::Matrix{Float64},
    linear_solver::String = "KLU",
)
    if linear_solver == "KLU"
        lodf_t = _calculate_LODF_matrix_KLU(a, ptdf)
    elseif linear_solver == "Dense"
        lodf_t = _calculate_LODF_matrix_DENSE(a, ptdf)
    elseif linear_solver == "MKLPardiso"
        if !USE_MKL
            error(
                "The MKL library is not available. Check that your hardware and operating system support MKL.",
            )
        end
        lodf_t = _calculate_LODF_matrix_MKLPardiso(a, ptdf)
    end
    return lodf_t
end

function _buildlodf(
    a::SparseArrays.SparseMatrixCSC{Int8, Int64},
    k::KLU.KLUFactorization{Float64, Int64},
    ba::SparseArrays.SparseMatrixCSC{Float64, Int64},
    ref_bus_positions::Set{Int},
    linear_solver::String,
)
    if linear_solver == "KLU"
        lodf_t = _calculate_LODF_matrix_KLU(a, k, ba, ref_bus_positions)
    else
        error("Other linear solvers are not implemented.")
    end
    return lodf_t
end

function _calculate_LODF_matrix_KLU(
    a::SparseArrays.SparseMatrixCSC{Int8, Int64},
    k::KLU.KLUFactorization{Float64, Int64},
    ba::SparseArrays.SparseMatrixCSC{Float64, Int64},
    ref_bus_positions::Set{Int64},
)
    linecount = size(ba, 2)
    # get inverse of aba
    first_ = zeros(size(a, 2), size(a, 1))
    valid_ix = setdiff(1:size(a, 2), ref_bus_positions)
    copyto!(first_, transpose(a))
    first_[valid_ix, :] = KLU.solve!(k, first_[valid_ix, :])
    first_[collect(ref_bus_positions), :] .= 0.0
    ptdf_denominator = first_' * ba

    m_I = Int[]
    m_V = Float64[]
    for iline in 1:linecount
        if (1.0 - ptdf_denominator[iline, iline]) < 1.0E-06
            push!(m_I, iline)
            push!(m_V, 1.0)
        else
            push!(m_I, iline)
            push!(m_V, 1 - ptdf_denominator[iline, iline])
        end
    end
    Dem_LU = klu(SparseArrays.sparse(m_I, m_I, m_V))
    KLU.solve!(Dem_LU, ptdf_denominator)
    ptdf_denominator[SparseArrays.diagind(ptdf_denominator)] .= -1.0
    return ptdf_denominator
end

function _calculate_LODF_matrix_KLU(
    a::SparseArrays.SparseMatrixCSC{Int8, Int},
    ptdf::Matrix{Float64},
)
    linecount = size(ptdf, 2)
    ptdf_denominator_t = a * ptdf
    m_I = Int[]
    m_V = Float64[]
    for iline in 1:linecount
        if (1.0 - ptdf_denominator_t[iline, iline]) < 1.0E-06
            push!(m_I, iline)
            push!(m_V, 1.0)
        else
            push!(m_I, iline)
            push!(m_V, 1 - ptdf_denominator_t[iline, iline])
        end
    end
    Dem_LU = klu(SparseArrays.sparse(m_I, m_I, m_V))
    lodf_t = Dem_LU \ ptdf_denominator_t
    lodf_t[SparseArrays.diagind(lodf_t)] .= -1.0

    return lodf_t
end

function _calculate_LODF_matrix_DENSE(
    a::SparseArrays.SparseMatrixCSC{Int8, Int},
    ptdf::Matrix{Float64},
)
    linecount = size(ptdf, 2)
    ptdf_denominator_t = a * ptdf
    m_V = Float64[]
    for iline in 1:linecount
        if (1.0 - ptdf_denominator_t[iline, iline]) < 1.0E-06
            push!(m_V, 1.0)
        else
            push!(m_V, 1.0 - ptdf_denominator_t[iline, iline])
        end
    end
    (mV, bipiv, binfo) = getrf!(Matrix(LinearAlgebra.diagm(m_V)))
    _binfo_check(binfo)
    getrs!('N', mV, bipiv, ptdf_denominator_t)
    ptdf_denominator_t[LinearAlgebra.diagind(ptdf_denominator_t)] .= -1.0
    return ptdf_denominator_t
end

function _pardiso_sequential_LODF!(
    lodf_t::Matrix{Float64},
    A::SparseArrays.SparseMatrixCSC{Float64, Int},
    ptdf_denominator_t::Matrix{Float64},
    chunk_size::Int = DEFAULT_LODF_CHUNK_SIZE,
)
    @info "Line Count too large for single compute using Pardiso. Employing Sequential Calulations using a chunk_size=$(chunk_size)"
    linecount = size(lodf_t, 1)
    @assert LinearAlgebra.ishermitian(A)
    ps = Pardiso.MKLPardisoSolver()
    Pardiso.set_matrixtype!(ps, Pardiso.REAL_SYM)
    Pardiso.pardisoinit(ps)
    # Pardiso.set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)
    defaults = Pardiso.get_iparms(ps)
    Pardiso.set_iparm!(ps, 1, 1)
    for (ix, v) in enumerate(defaults[2:end])
        Pardiso.set_iparm!(ps, ix + 1, v)
    end
    Pardiso.set_iparm!(ps, 2, 2)
    Pardiso.set_iparm!(ps, 59, 2)
    Pardiso.set_iparm!(ps, 12, 1)
    Pardiso.set_iparm!(ps, 11, 0)
    Pardiso.set_iparm!(ps, 13, 0)
    Pardiso.set_iparm!(ps, 32, 1)
    #Pardiso.set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)
    Pardiso.set_phase!(ps, Pardiso.ANALYSIS)
    Pardiso.pardiso(
        ps,
        lodf_t,
        A,
        ptdf_denominator_t,
    )

    Pardiso.set_phase!(ps, Pardiso.NUM_FACT)
    Pardiso.pardiso(
        ps,
        A,
        Float64[],
    )
    Pardiso.set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)
    i_count = 1
    tmp = zeros(Float64, linecount, chunk_size)
    while i_count <= linecount
        edge = min(i_count + chunk_size - 1, linecount)
        if linecount - edge <= 0
            tmp = tmp[:, 1:(edge - i_count + 1)]
        end
        Pardiso.pardiso(
            ps,
            tmp,
            A,
            ptdf_denominator_t[:, i_count:edge],
        )
        lodf_t[:, i_count:edge] .= tmp
        i_count = edge + 1
    end
    Pardiso.set_phase!(ps, Pardiso.RELEASE_ALL)
    Pardiso.pardiso(ps)
    return
end

function _pardiso_single_LODF!(
    lodf_t::Matrix{Float64},
    A::SparseArrays.SparseMatrixCSC{Float64, Int},
    ptdf_denominator_t::Matrix{Float64},
)
    @assert LinearAlgebra.ishermitian(A)
    ps = Pardiso.MKLPardisoSolver()
    Pardiso.set_matrixtype!(ps, Pardiso.REAL_SYM_POSDEF)
    Pardiso.pardisoinit(ps)
    Pardiso.set_iparm!(ps, 1, 1)
    defaults = Pardiso.get_iparms(ps)
    for (ix, v) in enumerate(defaults[2:end])
        Pardiso.set_iparm!(ps, ix + 1, v)
    end
    Pardiso.set_iparm!(ps, 2, 2)
    Pardiso.set_iparm!(ps, 59, 2)
    Pardiso.set_iparm!(ps, 12, 1)
    #Pardiso.set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)
    Pardiso.pardiso(
        ps,
        lodf_t,
        A,
        ptdf_denominator_t,
    )
    Pardiso.set_phase!(ps, Pardiso.RELEASE_ALL)
    Pardiso.pardiso(ps)
    return
end

function _calculate_LODF_matrix_MKLPardiso(
    a::SparseArrays.SparseMatrixCSC{Int8, Int},
    ptdf::Matrix{Float64},
)
    linecount = size(ptdf, 2)
    ptdf_denominator_t = a * ptdf
    m_I = Int[]
    m_V = Float64[]
    for iline in 1:linecount
        if (1.0 - ptdf_denominator_t[iline, iline]) < 1.0E-06
            push!(m_I, iline)
            push!(m_V, 1.0)
        else
            push!(m_I, iline)
            push!(m_V, 1 - ptdf_denominator_t[iline, iline])
        end
    end
    lodf_t = zeros(linecount, linecount)
    A = SparseArrays.sparse(m_I, m_I, m_V)
    if linecount > DEFAULT_LODF_CHUNK_SIZE
        _pardiso_sequential_LODF!(lodf_t, A, ptdf_denominator_t)
    else
        _pardiso_single_LODF!(lodf_t, A, ptdf_denominator_t)
    end
    lodf_t[LinearAlgebra.diagind(lodf_t)] .= -1.0
    return lodf_t
end

"""
Builds the LODF matrix given the data of branches and buses of the system.

# Arguments
- `branches`:
        vector of the System AC branches
- `buses::Vector{PSY.ACBus}`:
        vector of the System buses

# Keyword Arguments
- `linear_solver`::String
        linear solver to use for matrix evaluation.
        Available options: "KLU", "Dense" (OpenBLAS) or "MKLPardiso".
        Default solver: "KLU".
- `tol::Float64`:
        Tolerance to eliminate entries in the LODF matrix (default eps())
- `network_reduction::NetworkReduction`:
        Structure containing the details of the network reduction applied when computing the matrix
"""
function LODF(
    branches,
    buses::Vector{PSY.ACBus};
    linear_solver::String = "KLU",
    tol::Float64 = eps(),
    network_reduction::NetworkReduction = NetworkReduction(),
)
    # get axis names
    line_ax = [branch.name for branch in branches]
    axes = (line_ax, line_ax)
    line_map = make_ax_ref(line_ax)
    look_up = (line_map, line_map)
    bus_ax = [PSY.get_number(bus) for bus in buses]
    bus_lookup = make_ax_ref(bus_ax)
    # get network matrices
    ptdf_t, a = calculate_PTDF_matrix_KLU(branches, buses, bus_lookup, Float64[])

    if tol > eps()
        lodf_t = _buildlodf(a, ptdf_t, linear_solver)
        return LODF(
            sparsify(lodf_t, tol),
            axes,
            look_up,
            Ref(tol),
            network_reduction,
        )
    else
        return LODF(
            _buildlodf(a, ptdf_t, linear_solver),
            axes,
            look_up,
            Ref(tol),
            network_reduction,
        )
    end
end

"""
Builds the LODF matrix from a system.
Note that `network_reduction` kwargs is explicitly mentioned because needed inside of the function.

# Arguments
- `sys::PSY.System`:
        Power Systems system

# Keyword Arguments
- `network_reduction::NetworkReduction=NetworkReduction()`:
        Structure containing the details of the network reduction applied when computing the matrix
"""
function LODF(
    sys::PSY.System;
    network_reduction::NetworkReduction = NetworkReduction(),
    kwargs...,
)
    branches = get_ac_branches(sys, network_reduction.removed_branches)
    buses = get_buses(sys, network_reduction.bus_reduction_map)
    return LODF(branches, buses; network_reduction = network_reduction, kwargs...)
end

"""
Builds the LODF matrix given the Incidence Matrix and the PTDF matrix of the system.

NOTE: tol is referred to the LODF sparsification, not the PTDF one. PTDF matrix
must be considered as NON sparsified ("tol" argument not specified when calling
the PTDF method).

# Arguments
- `A::IncidenceMatrix`:
        Structure containing the Incidence matrix of the system.
- `PTDFm::PTDF`:
        Strucutre containing the transposed PTDF matrix of the system.
- `linear_solver::String`:
        Linear solver to be used. Options are "Dense" and "KLU".
- `tol::Float64`:
        Tolerance to eliminate entries in the LODF matrix (default eps()).
- `network_reduction::NetworkReduction`:
        Structure containing the details of the network reduction applied when computing the matrix
"""
function LODF(
    A::IncidenceMatrix,
    PTDFm::PTDF;
    linear_solver::String = "KLU",
    tol::Float64 = eps(),
    network_reduction::NetworkReduction = NetworkReduction(),
)
    validate_linear_solver(linear_solver)

    if PTDFm.tol.x > 1e-15
        warn_msg = string(
            "The argument `tol` in the PTDF matrix was set to a value different than the default one.\n",
            "The resulting LODF can include unexpected rounding errors.\n",
        )
        @warn(warn_msg)
        PTDFm_data = Matrix(PTDFm.data)
    else
        PTDFm_data = PTDFm.data
    end

    if !isequal(A.network_reduction, PTDFm.network_reduction)
        error("A and PTDF matrices have non-equivalent network reductions.")
    end
    ax_ref = make_ax_ref(A.axes[1])

    if tol > eps()
        lodf_t = _buildlodf(A.data, PTDFm_data, linear_solver)
        return LODF(
            sparsify(lodf_t, tol),
            (A.axes[1], A.axes[1]),
            (ax_ref, ax_ref),
            Ref(tol),
            network_reduction,
        )
    end
    return LODF(
        _buildlodf(A.data, PTDFm_data, linear_solver),
        (A.axes[1], A.axes[1]),
        (ax_ref, ax_ref),
        Ref(tol),
        network_reduction,
    )
end

"""
Builds the LODF matrix given the Incidence Matrix and the PTDF matrix of the system.

NOTE: this method does not support distributed slack bus.

# Arguments
- `A::IncidenceMatrix`:
        Structure containing the Incidence matrix of the system.
- `ABA::ABA_Matrix`:
        Structure containing the ABA matrix of the system.
- `BA::BA_Matrix`:
        Structure containing the transposed BA matrix of the system.
- `linear_solver::String`:
        Linear solver to be used. Options are "Dense" and "KLU".
- `tol::Float64`:
        Tolerance to eliminate entries in the LODF matrix (default eps()).
- `network_reduction::NetworkReduction`:
        Structure containing the details of the network reduction applied when computing the matrix
"""
function LODF(
    A::IncidenceMatrix,
    ABA::ABA_Matrix,
    BA::BA_Matrix;
    linear_solver::String = "KLU",
    tol::Float64 = eps(),
    network_reduction::NetworkReduction = NetworkReduction(), #TODO The public API should not have a way to pass in network_reduction to this function which uses existing matrices to compute the LODF. 
)
    if !(
        isequal(A.network_reduction, BA.network_reduction) &&
        isequal(BA.network_reduction, ABA.network_reduction)
    )
        error(
            "Mismatch in `NetworkReduction`, A, BA, and ABA matrices must be computed with the same network reduction.",
        )
    end
    validate_linear_solver(linear_solver)
    ax_ref = make_ax_ref(A.axes[1])

    if tol > eps()
        lodf_t = _buildlodf(A.data, ABA.K, BA.data, A.ref_bus_positions, linear_solver)
        return LODF(
            sparsify(lodf_t, tol),
            (A.axes[1], A.axes[1]),
            (ax_ref, ax_ref),
            Ref(tol),
            network_reduction,
        )
    end
    return LODF(
        _buildlodf(A.data, ABA.K, BA.data, A.ref_bus_positions, linear_solver),
        (A.axes[1], A.axes[1]),
        (ax_ref, ax_ref),
        Ref(tol),
        network_reduction,
    )
end

############################################################
# auxiliary functions for getting data from LODF structure #
############################################################

# NOTE: the LODF matrix is saved as transposed!

function Base.getindex(A::LODF, selected_line, outage_line)
    i, j = to_index(A, outage_line, selected_line)
    return A.data[i, j]
end

function Base.getindex(
    A::LODF,
    selected_line_number::Union{Int, Colon},
    outage_line_number::Union{Int, Colon},
)
    return A.data[outage_line_number, selected_line_number]
end

function get_lodf_data(lodf::LODF)
    return transpose(lodf.data)
end

function get_branch_ax(lodf::LODF)
    return lodf.axes[1]
end

function get_tol(lodf::LODF)
    return lodf.tol
end
