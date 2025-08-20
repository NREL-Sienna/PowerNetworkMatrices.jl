"""
Nodal admittance matrix (Ybus) is an N x N matrix describing a power system with N buses. It represents the nodal admittance of the buses in a power system.

The Ybus Struct is indexed using the Bus Numbers, no need for them to be sequential

The fields yft and ytf are the branch admittance matrices for the from-to and to-from branch admittances respectively. The rows correspond to branches and the columns to buses.
The matrix columns are mapped to buses using fb, tb arrays of the matrix columns that correspond to the `from` and `to` buses.
Using yft, ytf, and the voltage vector, the branch currents and power flows can be calculated.
"""
struct Ybus{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{ComplexF32}
    data::SparseArrays.SparseMatrixCSC{ComplexF32, Int}
    adjacency_data::SparseArrays.SparseMatrixCSC{Int8, Int}
    axes::Ax
    lookup::L
    subnetwork_axes::Dict{Int, Ax}
    arc_subnetwork_axis::Dict{Int, Vector{Tuple{Int, Int}}}
    network_reduction_data::NetworkReductionData
    branch_admittance_from_to::Union{ArcAdmittanceMatrix, Nothing}
    branch_admittance_to_from::Union{ArcAdmittanceMatrix, Nothing}
end

get_axes(M::Ybus) = M.axes
get_lookup(M::Ybus) = M.lookup
get_ref_bus(M::Ybus) = sort!(collect(keys(M.subnetwork_axes)))
get_ref_bus_position(M::Ybus) = [get_bus_lookup(M)[x] for x in keys(M.subnetwork_axes)]
get_network_reduction_data(M::Ybus) = M.network_reduction_data
get_bus_axis(M::Ybus) = M.axes[1]
get_bus_lookup(M::Ybus) = M.lookup[1]

function get_isolated_buses(M::Ybus)
    return [x for x in keys(M.subnetwork_axes) if length(M.subnetwork_axes[x][1]) == 1]
end

function get_default_reduction(sys::PSY.System)
    ybus = Ybus(sys)
    return ybus.network_reduction_data
end

function get_reduction(
    ybus::Ybus,
    sys::PSY.System,
    reduction::RadialReduction,
)
    A = IncidenceMatrix(ybus)
    return get_reduction(A, sys, reduction)
end

function add_to_branch_maps!(nr::NetworkReductionData, arc::PSY.Arc, br::PSY.Branch)
    direct_branch_map = get_direct_branch_map(nr)
    reverse_direct_branch_map = get_reverse_direct_branch_map(nr)
    parallel_branch_map = get_parallel_branch_map(nr)
    reverse_parallel_branch_map = get_reverse_parallel_branch_map(nr)
    arc_tuple = get_arc_tuple(arc, nr)
    if haskey(parallel_branch_map, arc_tuple)
        push!(parallel_branch_map[arc_tuple], br)
        reverse_parallel_branch_map[br] = arc_tuple
    elseif haskey(direct_branch_map, arc_tuple)
        corresponding_branch = direct_branch_map[arc_tuple]
        delete!(direct_branch_map, arc_tuple)
        delete!(reverse_direct_branch_map, corresponding_branch)
        parallel_branch_map[arc_tuple] = Set([corresponding_branch, br])
        reverse_parallel_branch_map[corresponding_branch] = arc_tuple
        reverse_parallel_branch_map[br] = arc_tuple
    else
        direct_branch_map[arc_tuple] = br
        reverse_direct_branch_map[br] = arc_tuple
    end
    return
end

function add_to_branch_maps!(
    nr::NetworkReductionData,
    primary_star_arc::PSY.Arc,
    secondary_star_arc::PSY.Arc,
    tertiary_star_arc::PSY.Arc,
    br::PSY.ThreeWindingTransformer,
)
    transformer3W_map = get_transformer3W_map(nr)
    reverse_transformer3W_map = get_reverse_transformer3W_map(nr)
    if PSY.get_available_primary(br)
        primary_star_arc_tuple = get_arc_tuple(primary_star_arc, nr)
        transformer3W_map[primary_star_arc_tuple] = (br, 1)
        reverse_transformer3W_map[(br, 1)] = primary_star_arc_tuple
    end
    if PSY.get_available_secondary(br)
        secondary_star_arc_tuple = get_arc_tuple(secondary_star_arc, nr)
        transformer3W_map[secondary_star_arc_tuple] = (br, 2)
        reverse_transformer3W_map[(br, 2)] = secondary_star_arc_tuple
    end
    if PSY.get_available_tertiary(br)
        tertiary_star_arc_tuple = get_arc_tuple(tertiary_star_arc, nr)
        transformer3W_map[tertiary_star_arc_tuple] = (br, 3)
        reverse_transformer3W_map[(br, 3)] = tertiary_star_arc_tuple
    end
    return
end

function add_branch_entries_to_ybus!(
    y11::Vector{ComplexF32},
    y12::Vector{ComplexF32},
    y21::Vector{ComplexF32},
    y22::Vector{ComplexF32},
    branch_ix::Int,
    br::PSY.ACTransmission,
)
    Y11, Y12, Y21, Y22 = ybus_branch_entries(br)
    y11[branch_ix] = Y11
    y12[branch_ix] = Y12
    y21[branch_ix] = Y21
    y22[branch_ix] = Y22
    return
end

function add_branch_entries_to_indexing_maps!(
    num_bus::Dict{Int, Int},
    branch_ix::Int,
    nr::NetworkReductionData,
    fb::Vector{Int},
    tb::Vector{Int},
    adj::SparseArrays.SparseMatrixCSC{Int8, Int},
    br::PSY.ACTransmission,
)
    arc = PSY.get_arc(br)
    add_to_branch_maps!(nr, arc, br)
    bus_from_no, bus_to_no = get_bus_indices(arc, num_bus, nr)
    adj[bus_from_no, bus_to_no] = 1
    adj[bus_to_no, bus_from_no] = -1
    fb[branch_ix] = bus_from_no
    tb[branch_ix] = bus_to_no
    return
end

_get_shunt(br::PSY.ACTransmission, node::Symbol) =
    PSY.get_g(br)[node] + 1im * PSY.get_b(br)[node]
_get_shunt(::PSY.DiscreteControlledACBranch, ::Symbol) = zero(ComplexF32)

"""Ybus entries for a `Line` or a `DiscreteControlledACBranch`."""
function ybus_branch_entries(br::PSY.ACTransmission)
    Y_l = (1 / (PSY.get_r(br) + PSY.get_x(br) * 1im))
    Y11 = Y_l + _get_shunt(br, :from)
    if !isfinite(Y11) || !isfinite(Y_l)
        error(
            "Data in $(PSY.get_name(br)) is incorrect. r = $(PSY.get_r(br)), x = $(PSY.get_x(br))",
        )
    end
    Y12 = -Y_l
    Y21 = Y12
    Y22 = Y_l + _get_shunt(br, :to)
    return (Y11, Y12, Y21, Y22)
end

_get_tap(::PSY.Transformer2W) = one(ComplexF32)
_get_tap(br::PSY.TwoWindingTransformer) = PSY.get_tap(br)

"""Ybus entries for a `Transformer2W`, `TapTransformer`, or `PhaseShiftingTransformer`."""
function ybus_branch_entries(br::PSY.TwoWindingTransformer)
    Y_t = 1 / (PSY.get_r(br) + PSY.get_x(br) * 1im)
    tap = _get_tap(br) * exp(PSY.get_α(br) * 1im)
    c_tap = _get_tap(br) * exp(-1 * PSY.get_α(br) * 1im)
    y_shunt = PSY.get_primary_shunt(br)
    Y11 = Y_t / abs2(tap)
    if !isfinite(Y11) || !isfinite(Y_t) || !isfinite(y_shunt * c_tap)
        error(
            "Data in $(summary(br)) gives a non-finite Ybus entry; check input data.",
        )
    end
    Y12 = -Y_t / c_tap
    Y21 = -Y_t / tap
    Y22 = Y_t
    return (Y11 + y_shunt, Y12, Y21, Y22)
end

"""Ybus branch entries for an arc in the wye model of a `ThreeWindingTransformer`."""
function ybus_branch_entries(tp::Tuple{PSY.ThreeWindingTransformer, Int})
    (br, winding_number) = tp
    if winding_number == 1
        Y_t = 1 / (PSY.get_r_primary(br) + PSY.get_x_primary(br) * 1im)
        tap = PSY.get_primary_turns_ratio(br) * exp(PSY.get_α_primary(br) * 1im)
        c_tap = PSY.get_primary_turns_ratio(br) * exp(-1 * PSY.get_α_primary(br) * 1im)
        Y11 = Y_t / abs2(tap)
        y_shunt = PSY.get_g(br) + im * PSY.get_b(br)
        if !isfinite(Y11) || !isfinite(y_shunt)
            error(
                "Data in $(PSY.get_name(br)) is incorrect.
                r_p = $(PSY.get_r_primary(br)), x_p = $(PSY.get_x_primary(br)), primary_turns_ratio = $(PSY.get_primary_turns_ratio(br))",
            )
        end
        # primary bus alone gets the shunt term
        Y11 += y_shunt
    elseif winding_number == 2
        Y_t = 1 / (PSY.get_r_secondary(br) + PSY.get_x_secondary(br) * 1im)
        tap = PSY.get_secondary_turns_ratio(br) * exp(PSY.get_α_secondary(br) * 1im)
        c_tap = PSY.get_secondary_turns_ratio(br) * exp(-1 * PSY.get_α_secondary(br) * 1im)
        Y11 = Y_t / abs2(tap)
        if !isfinite(Y11)
            error(
                "Data in $(PSY.get_name(br)) is incorrect.
                r_s = $(PSY.get_r_secondary(br)), x_s = $(PSY.get_x_secondary(br)), secondary_turns_ratio = $(PSY.get_secondary_turns_ratio(br))",
            )
        end
    elseif winding_number == 3
        Y_t = 1 / (PSY.get_r_tertiary(br) + PSY.get_x_tertiary(br) * 1im)
        tap = PSY.get_tertiary_turns_ratio(br) * exp(PSY.get_α_tertiary(br) * 1im)
        c_tap = PSY.get_tertiary_turns_ratio(br) * exp(-1 * PSY.get_α_tertiary(br) * 1im)
        Y11 = Y_t / abs2(tap)
        if !isfinite(Y11)
            error(
                "Data in $(PSY.get_name(br)) is incorrect.
                r_t = $(PSY.get_r_tertiary(br)), x_t = $(PSY.get_x_tertiary(br)), tertiary_turns_ratio = $(PSY.get_tertiary_turns_ratio(br))",
            )
        end
    end
    Y12 = (-Y_t / c_tap)
    Y21 = (-Y_t / tap)
    Y22 = Y_t
    return (Y11, Y12, Y21, Y22)
end

"""Handles ybus entries for most 2-node AC branches. The types handled here are:
`Line`, `DiscreteControlledACBranch`, `Transformer2W`, `TapTransformer`, and `PhaseShiftingTransformer`.
"""
function _ybus!(
    y11::Vector{ComplexF32},
    y12::Vector{ComplexF32},
    y21::Vector{ComplexF32},
    y22::Vector{ComplexF32},
    br::PSY.ACTransmission,
    num_bus::Dict{Int, Int},
    branch_ix::Int,
    fb::Vector{Int},
    tb::Vector{Int},
    nr::NetworkReductionData,
    adj::SparseArrays.SparseMatrixCSC{Int8, Int},
)
    add_branch_entries_to_indexing_maps!(num_bus, branch_ix, nr, fb, tb, adj, br)
    add_branch_entries_to_ybus!(y11, y12, y21, y22, branch_ix, br)
    return
end

function _ybus!(
    y11::Vector{ComplexF32},
    y12::Vector{ComplexF32},
    y21::Vector{ComplexF32},
    y22::Vector{ComplexF32},
    br::PSY.DynamicBranch,
    num_bus::Dict{Int, Int},
    branch_ix::Int,
    fb::Vector{Int},
    tb::Vector{Int},
    nr::NetworkReductionData,
    adj::SparseArrays.SparseMatrixCSC{Int8, Int},
)
    _ybus!(
        y11,
        y12,
        y21,
        y22,
        br.branch,
        num_bus,
        branch_ix,
        fb,
        tb,
        nr,
        adj,
    )
    return
end

function _ybus!(
    y11::Vector{ComplexF32},
    y12::Vector{ComplexF32},
    y21::Vector{ComplexF32},
    y22::Vector{ComplexF32},
    br::PSY.ThreeWindingTransformer,
    num_bus::Dict{Int, Int},
    offset_ix::Int,
    fb::Vector{Int},
    tb::Vector{Int},
    ix::Int,
    nr::NetworkReductionData,
    adj::SparseArrays.SparseMatrixCSC{Int8, Int},
)
    primary_star_arc = PSY.get_primary_star_arc(br)
    secondary_star_arc = PSY.get_secondary_star_arc(br)
    tertiary_star_arc = PSY.get_tertiary_star_arc(br)
    add_to_branch_maps!(nr, primary_star_arc, secondary_star_arc, tertiary_star_arc, br)    #TODO - check this for case of some arcs unavailable
    primary_available = PSY.get_available_primary(br)
    secondary_available = PSY.get_available_secondary(br)
    tertiary_available = PSY.get_available_tertiary(br)
    n_entries = 0
    if primary_available
        primary_ix, star_ix = get_bus_indices(primary_star_arc, num_bus, nr)
        adj[primary_ix, star_ix] = 1
        adj[star_ix, primary_ix] = -1
        fb[offset_ix + ix + n_entries] = primary_ix
        tb[offset_ix + ix + n_entries] = star_ix
        (Y11, Y12, Y21, Y22) = ybus_branch_entries((br, 1))
        y11[offset_ix + ix + n_entries] = Y11
        y12[offset_ix + ix + n_entries] = Y12
        y21[offset_ix + ix + n_entries] = Y21
        y22[offset_ix + ix + n_entries] = Y22
        n_entries += 1
    end
    if secondary_available
        secondary_ix, star_ix = get_bus_indices(secondary_star_arc, num_bus, nr)
        adj[secondary_ix, star_ix] = 1
        adj[star_ix, secondary_ix] = -1
        fb[offset_ix + ix + n_entries] = secondary_ix
        tb[offset_ix + ix + n_entries] = star_ix
        (Y11, Y12, Y21, Y22) = ybus_branch_entries((br, 2))
        y11[offset_ix + ix + n_entries] = Y11
        y12[offset_ix + ix + n_entries] = Y12
        y21[offset_ix + ix + n_entries] = Y21
        y22[offset_ix + ix + n_entries] = Y22
        n_entries += 1
    end
    if tertiary_available
        tertiary_ix, star_ix = get_bus_indices(tertiary_star_arc, num_bus, nr)
        adj[tertiary_ix, star_ix] = 1
        adj[star_ix, tertiary_ix] = -1
        fb[offset_ix + ix + n_entries] = tertiary_ix
        tb[offset_ix + ix + n_entries] = star_ix
        (Y11, Y12, Y21, Y22) = ybus_branch_entries((br, 3))
        y11[offset_ix + ix + n_entries] = Y11
        y12[offset_ix + ix + n_entries] = Y12
        y21[offset_ix + ix + n_entries] = Y21
        y22[offset_ix + ix + n_entries] = Y22
        n_entries += 1
    end
    return n_entries
end

function _ybus!(
    ysh::Vector{ComplexF32},
    fa::PSY.FixedAdmittance,
    num_bus::Dict{Int, Int},
    fa_ix::Int,
    sb::Vector{Int},
    nr::NetworkReductionData,
)
    bus_no = get_bus_index(fa, num_bus, nr)
    Y = PSY.get_Y(fa)
    sb[fa_ix] = bus_no
    if !isfinite(Y)
        error(
            "Data in $(PSY.get_name(fa)) is incorrect. Y = $(Y)",
        )
    end
    ysh[fa_ix] = Y
    return
end

#Note - PSSE does not include switched admittances in ymatrix
function _ybus!(
    ysh::Vector{ComplexF32},
    fa::PSY.SwitchedAdmittance,
    num_bus::Dict{Int, Int},
    fa_ix::Int,
    sb::Vector{Int},
    nr::NetworkReductionData,
)
    bus_no = get_bus_index(fa, num_bus, nr)
    sb[fa_ix] = bus_no
    ysh[fa_ix] = 0.0
    return
end

function _ybus!(
    ysh::Vector{ComplexF32},
    fa::PSY.StandardLoad,
    num_bus::Dict{Int, Int},
    fa_ix::Int,
    sb::Vector{Int},
    nr::NetworkReductionData,
)
    bus_no = get_bus_index(fa, num_bus, nr)
    Y = PSY.get_impedance_active_power(fa) - im * PSY.get_impedance_reactive_power(fa)
    if !isfinite(Y)
        error(
            "Data in $(PSY.get_name(fa)) is incorrect. Y = $(Y)",
        )
    end
    sb[fa_ix] = bus_no
    ysh[fa_ix] = Y
    return
end

function _buildybus!(
    network_reduction_data::NetworkReductionData,
    adj::SparseArrays.SparseMatrixCSC{Int8, Int},
    branches,
    transformer_3w::Vector{PSY.ThreeWindingTransformer},
    num_bus::Dict{Int, Int},
    fixed_admittances::Vector{PSY.FixedAdmittance},
    switched_admittances::Vector{PSY.SwitchedAdmittance},
    standard_loads::Vector{PSY.StandardLoad},
)
    branch_entries_transformer_3w = 0
    for br in transformer_3w
        branch_entries_transformer_3w += count([
            PSY.get_available_primary(br),
            PSY.get_available_secondary(br),
            PSY.get_available_tertiary(br),
        ])
    end
    branchcount = length(branches) + branch_entries_transformer_3w
    branchcount_no_3w = length(branches)
    fa_count = length(fixed_admittances)
    sa_count = length(switched_admittances)
    sl_count = length(standard_loads)
    fb = zeros(Int, branchcount)
    tb = zeros(Int, branchcount)
    sb = zeros(Int, fa_count + sa_count + sl_count)

    y11 = zeros(ComplexF32, branchcount)
    y12 = zeros(ComplexF32, branchcount)
    y21 = zeros(ComplexF32, branchcount)
    y22 = zeros(ComplexF32, branchcount)
    ysh = zeros(ComplexF32, fa_count + sa_count + sl_count)

    for (ix, b) in enumerate(branches)
        if PSY.get_name(b) == "init"
            throw(DataFormatError("The data in Branch is invalid"))
        end
        _ybus!(y11, y12, y21, y22, b, num_bus, ix, fb, tb, network_reduction_data, adj)
    end

    ix = 1
    for b in transformer_3w
        if PSY.get_name(b) == "init"
            throw(DataFormatError("The data in Transformer3W is invalid"))
        end
        n_entries = _ybus!(
            y11,
            y12,
            y21,
            y22,
            b,
            num_bus,
            branchcount_no_3w,
            fb,
            tb,
            ix,
            network_reduction_data,
            adj,
        )
        ix += n_entries
    end
    for (ix, fa) in enumerate([fixed_admittances; switched_admittances; standard_loads])
        _ybus!(ysh, fa, num_bus, ix, sb, network_reduction_data)
    end
    return (
        y11,
        y12,
        y21,
        y22,
        ysh,
        fb,
        tb,
        sb,
    )
end

"""
Builds a Ybus from the system. The return is a Ybus Array indexed with the bus numbers and the arc tuples.
"""
function Ybus(
    sys::PSY.System;
    make_branch_admittance_matrices::Bool = false,
    network_reductions::Vector{NetworkReduction} = NetworkReduction[],
    include_constant_impedance_loads = true,
    subnetwork_algorithm = iterative_union_find,
    kwargs...,
)
    ref_bus_numbers = Set{Int}()
    nr = NetworkReductionData()
    bus_reduction_map = get_bus_reduction_map(nr)
    reverse_bus_search_map = get_reverse_bus_search_map(nr)

    #Checking for isolated buses; building bus map.
    for b in PSY.get_components(x -> PSY.get_available(x), PSY.ACBus, sys)
        if PSY.get_bustype(b) != ACBusTypes.ISOLATED
            bus_reduction_map[PSY.get_number(b)] = Set{Int}()
            if PSY.get_bustype(b) == ACBusTypes.REF
                push!(ref_bus_numbers, PSY.get_number(b))
            end
        else
            @debug "Found available isolated bus $(PSY.get_name(b)) with number $(PSY.get_number(b)). This is excluded from the Ybus build."
            push!(nr.removed_buses, PSY.get_number(b))
        end
    end

    #Building map for removed Breaker/Switches
    breaker_switches = Vector{PSY.DiscreteControlledACBranch}()
    for br in
        PSY.get_components(x -> PSY.get_available(x), PSY.DiscreteControlledACBranch, sys)
        r = PSY.get_r(br)
        x = PSY.get_x(br)
        status = PSY.get_branch_status(br)
        if status == PSY.DiscreteControlledBranchStatus.CLOSED
            if r == 0.0 && x < ZERO_IMPEDANCE_LINE_REACTANCE_THRESHOLD
                from_bus_number = PSY.get_number(PSY.get_from(PSY.get_arc(br)))
                to_bus_number = PSY.get_number(PSY.get_to(PSY.get_arc(br)))
                if haskey(reverse_bus_search_map, from_bus_number)
                    reduced_from_bus_number = reverse_bus_search_map[from_bus_number]
                    push!(
                        get!(bus_reduction_map, reduced_from_bus_number, Set{Int}),
                        to_bus_number,
                    )
                    delete!(bus_reduction_map, to_bus_number)
                    reverse_bus_search_map[to_bus_number] = reduced_from_bus_number
                else
                    s1 = get(bus_reduction_map, from_bus_number, Set{Int}())
                    s2 = union(
                        get(bus_reduction_map, to_bus_number, Set{Int}(to_bus_number)),
                        to_bus_number,
                    )
                    bus_reduction_map[from_bus_number] = union(s1, s2)
                    delete!(bus_reduction_map, to_bus_number)
                    for x in s2
                        reverse_bus_search_map[x] = from_bus_number
                    end
                end
                push!(nr.removed_arcs, (from_bus_number, to_bus_number))
            else
                push!(breaker_switches, br)
            end
        end
    end

    bus_ax = sort!([x for x in keys(bus_reduction_map)])
    axes = (bus_ax, bus_ax)
    bus_lookup = Dict{Int, Int}()
    lookup = (bus_lookup, bus_lookup)
    busnumber = length(bus_ax)
    for (ix, b) in enumerate(bus_ax)
        bus_lookup[b] = ix
    end
    adj = SparseArrays.spdiagm(ones(Int8, busnumber))
    branches = collect(
        PSY.get_components(
            x ->
                PSY.get_available(x) &&
                    typeof(x) ∉ [
                        PSY.Transformer3W,
                        PSY.PhaseShiftingTransformer3W,
                        PSY.DiscreteControlledACBranch,
                    ],
            PSY.ACTransmission,
            sys,
        ),
    )
    branches = vcat(branches, breaker_switches)
    transformer_3W =
        collect(
            PSY.get_components(
                x -> PSY.get_available(x),
                PSY.ThreeWindingTransformer,
                sys,
            ),
        )
    fixed_admittances =
        collect(PSY.get_components(x -> PSY.get_available(x), PSY.FixedAdmittance, sys))
    switched_admittances =
        collect(PSY.get_components(x -> PSY.get_available(x), PSY.SwitchedAdmittance, sys))
    standard_loads = if include_constant_impedance_loads
        collect(PSY.get_components(x -> PSY.get_available(x), PSY.StandardLoad, sys))
    else
        PSY.StandardLoad[]
    end
    y11, y12, y21, y22, ysh, fb, tb, sb =
        _buildybus!(
            nr,
            adj,
            branches,
            transformer_3W,
            bus_lookup,
            fixed_admittances,
            switched_admittances,
            standard_loads,
        )
    ybus = SparseArrays.sparse(
        [fb; fb; tb; tb; sb],  # row indices
        [fb; tb; fb; tb; sb],  # column indices
        [y11; y12; y21; y22; ysh],  # values
        busnumber,  # size (rows) - setting this explicitly is necessary for the case there are no branches
        busnumber,  # size (columns) - setting this explicitly is necessary for the case there are no branches
    )
    SparseArrays.dropzeros!(ybus)

    if make_branch_admittance_matrices
        arc_axis = get_arc_axis(fb, tb, bus_ax)
        arc_lookup = Dict{Tuple{Int, Int}, Int}()
        for (ix, arc_tuple) in enumerate(arc_axis)
            arc_lookup[arc_tuple] = ix
        end
        rows_ix = [arc_lookup[(x, y)] for (x, y) in zip(bus_ax[fb], bus_ax[tb])]
        yft_data = SparseArrays.sparse(
            vcat(rows_ix, rows_ix),
            [fb; tb],
            [y11; y12],
            length(arc_axis),
            busnumber,
        )
        ytf_data = SparseArrays.sparse(
            vcat(rows_ix, rows_ix),
            [tb; fb],
            [y22; y21],
            length(arc_axis),
            busnumber,
        )
        branch_admittance_from_to = ArcAdmittanceMatrix(
            yft_data,
            (arc_axis, bus_ax),
            (arc_lookup, bus_lookup),
            nr,
            :FromTo,
        )
        branch_admittance_to_from = ArcAdmittanceMatrix(
            ytf_data,
            (arc_axis, bus_ax),
            (arc_lookup, bus_lookup),
            nr,
            :ToFrom,
        )
    else
        branch_admittance_from_to = nothing
        branch_admittance_to_from = nothing
    end
    if length(bus_lookup) > 1
        subnetworks = assign_reference_buses!(
            find_subnetworks(ybus, bus_ax; subnetwork_algorithm = subnetwork_algorithm),
            ref_bus_numbers,
        )
        if length(subnetworks) > 1
            @warn "More than one island found; Network is not connected"
        end
    else
        subnetworks = Dict{Int, Set{Int}}(bus_lookup[1] => Set(bus_ax))
    end
    subnetwork_axes = _make_bus_subnetwork_axes(subnetworks)
    arc_subnetwork_axis = _make_arc_subnetwork_axis(subnetworks, nr)
    ybus = Ybus(
        ybus,
        adj,
        axes,
        lookup,
        subnetwork_axes,
        arc_subnetwork_axis,
        nr,
        branch_admittance_from_to,
        branch_admittance_to_from,
    )

    for nr in network_reductions
        ybus = build_reduced_ybus(ybus, sys, nr)
    end
    return ybus
end

#TODO - handle arc axis consistently between BranchAdmittanceMatrices and IncidenceMatrix
function get_arc_axis(fb::Vector{Int}, tb::Vector{Int}, bus_axis::Vector{Int})
    return unique(collect(zip(bus_axis[fb], bus_axis[tb])))
end

"""
Make subnetwork axes for BA_Matrix
"""
function make_bus_arc_subnetwork_axes(ybus::Ybus)
    subnetwork_axes = Dict{Int, Tuple{Vector{Int}, Vector{Tuple{Int, Int}}}}()
    for key in keys(ybus.subnetwork_axes)
        subnetwork_axes[key] = (ybus.subnetwork_axes[key][1], ybus.arc_subnetwork_axis[key])
    end
    return subnetwork_axes
end

"""
Make subnetwork axes for IncidenceMatrix
"""
function make_arc_bus_subnetwork_axes(ybus::Ybus)
    subnetwork_axes = Dict{Int, Tuple{Vector{Tuple{Int, Int}}, Vector{Int}}}()
    for key in keys(ybus.subnetwork_axes)
        subnetwork_axes[key] = (ybus.arc_subnetwork_axis[key], ybus.subnetwork_axes[key][1])
    end
    return subnetwork_axes
end

function _make_bus_subnetwork_axes(subnetworks::Dict{Int, Set{Int}})
    subnetwork_axes = Dict{Int, Tuple{Vector{Int}, Vector{Int}}}()
    for (k, v) in subnetworks
        subnetwork_axes[k] = (collect(v), collect(v))
    end
    return subnetwork_axes
end

function _make_arc_subnetwork_axis(
    subnetworks::Dict{Int, Set{Int}},
    nr::NetworkReductionData,
)
    arc_ax = get_arc_axis(nr)
    arc_subnetwork_axis = Dict{Int, Vector{Tuple{Int, Int}}}()
    for k in keys(subnetworks)
        arc_subnetwork_axis[k] = Vector{Tuple{Int, Int}}()
    end
    for arc in arc_ax
        for (k, v) in subnetworks
            if arc[1] ∈ v || arc[2] in v
                subnetwork = get!(arc_subnetwork_axis, k, Vector{Tuple{Int, Int}}())
                push!(subnetwork, arc)
                break
            end
        end
    end
    return arc_subnetwork_axis
end

function build_reduced_ybus(
    ybus::Ybus,
    sys::PSY.System,
    network_reduction::NetworkReduction,
)
    network_reduction_data = get_reduction(ybus, sys, network_reduction)
    return _apply_reduction(ybus, network_reduction_data)
end

#NOTE: this is the key function that composes sequential reductions; this function needs cleanup, review, and more testing.
function _apply_reduction(ybus::Ybus, nr_new::NetworkReductionData)
    #TODO - we should change this to do this validation before we compute the reduction
    validate_reduction_type(
        get_reductions(nr_new),
        get_reductions(get_network_reduction_data(ybus)),
    )
    remake_reverse_direct_branch_map = false
    remake_reverse_parallel_branch_map = false
    remake_reverse_series_branch_map = false
    remake_reverse_transformer3W_map = false
    data = get_data(ybus)
    adjacency_data = ybus.adjacency_data
    lookup = get_lookup(ybus)
    bus_lookup = lookup[1]
    nr = ybus.network_reduction_data
    bus_numbers_to_remove = Vector{Int}()
    for (k, v) in nr_new.reverse_bus_search_map
        nr.reverse_bus_search_map[k] = v
        push!(bus_numbers_to_remove, k)
    end
    for x in nr_new.removed_buses
        push!(nr.removed_buses, x)
        push!(bus_numbers_to_remove, x)
        delete!(nr.bus_reduction_map, x)
    end
    for (k, v) in nr_new.bus_reduction_map
        if haskey(nr.bus_reduction_map, k)
            union!(nr.bus_reduction_map[k], nr_new.bus_reduction_map[k])
            for x in v
                delete!(nr.bus_reduction_map, x)
            end
        else
            error("Bus $k was previously reduced")
        end
    end
    for x in nr_new.removed_arcs
        push!(nr.removed_arcs, x)
        if haskey(nr.direct_branch_map, x)
            remake_reverse_direct_branch_map = true
            pop!(nr.direct_branch_map, x)
        elseif haskey(nr.parallel_branch_map, x)
            remake_reverse_parallel_branch_map = true
            pop!(nr.parallel_branch_map, x)
        elseif haskey(nr.series_branch_map, x)
            remake_reverse_series_branch_map = true
            pop!(nr.series_branch_map, x)
        elseif haskey(nr.transformer3W_map, x)
            remake_reverse_transformer3W_map = true
            pop!(nr.transformer3W_map, x)
        end
    end
    # Add additional entries to the ybus corresponding to the equivalent series arcs
    new_y_ft, new_y_tf = _add_series_branches_to_ybus!(
        ybus.data,
        get_bus_lookup(ybus),
        ybus.branch_admittance_from_to,
        ybus.branch_admittance_to_from,
        nr_new.series_branch_map,
        nr,
    )

    remake_reverse_direct_branch_map && _remake_reverse_direct_branch_map!(nr)
    remake_reverse_parallel_branch_map && _remake_reverse_parallel_branch_map!(nr)
    remake_reverse_series_branch_map && _remake_reverse_series_branch_map!(nr)
    remake_reverse_transformer3W_map && _remake_reverse_transformer3W_map!(nr)
    #Assumes only the last reduction (Ward) can add branches and admittances
    nr.added_branch_map = nr_new.added_branch_map
    nr.added_admittance_map = nr_new.added_admittance_map
    if isempty(nr.series_branch_map)
        nr.series_branch_map = nr_new.series_branch_map
        nr.reverse_series_branch_map = nr_new.reverse_series_branch_map
    elseif !isempty(nr_new.series_branch_map) && !isempty(nr.series_branch_map)
        error(
            "Cannot compose series branch maps; should not apply multiple reductions that generate series branch maps",
        )
    end
    for (bus_no, admittance) in nr.added_admittance_map
        data[bus_lookup[bus_no], bus_lookup[bus_no]] += admittance
    end
    for (bus_tuple, admittance) in nr.added_branch_map
        bus_from, bus_to = bus_tuple
        data[bus_lookup[bus_from], bus_lookup[bus_to]] += admittance
        data[bus_lookup[bus_to], bus_lookup[bus_from]] += admittance
    end
    add_reduction!(nr.reductions, nr_new.reductions)
    union!(nr.irreducible_buses, nr_new.irreducible_buses)
    bus_ax = setdiff(get_bus_axis(ybus), bus_numbers_to_remove)
    bus_lookup = make_ax_ref(bus_ax)
    bus_ix = [ybus.lookup[1][x] for x in bus_ax]
    adjacency_data = adjacency_data[bus_ix, bus_ix]
    data = data[bus_ix, bus_ix]

    subnetwork_axes, arc_subnetwork_axis =
        _make_subnetwork_axes(ybus, bus_numbers_to_remove, nr_new.removed_arcs)

    if new_y_ft !== nothing
        arc_ax = setdiff(get_arc_axis(new_y_ft), nr_new.removed_arcs)
        arc_remove_ixs = indexin(nr_new.removed_arcs, get_arc_axis(new_y_ft))
        arc_keep_ixs = setdiff(collect(1:length(get_arc_axis(new_y_ft))), arc_remove_ixs)
        arc_lookup = make_ax_ref(arc_ax)
        yft_data = new_y_ft.data[arc_keep_ixs, bus_ix]
        ytf_data = new_y_tf.data[arc_keep_ixs, bus_ix]

        branch_admittance_from_to = ArcAdmittanceMatrix(
            yft_data,
            (arc_ax, bus_ax),
            (arc_lookup, bus_lookup),
            nr,
            :FromTo,
        )
        branch_admittance_to_from = ArcAdmittanceMatrix(
            ytf_data,
            (arc_ax, bus_ax),
            (arc_lookup, bus_lookup),
            nr,
            :ToFrom,
        )
    else
        branch_admittance_from_to = ybus.branch_admittance_from_to
        branch_admittance_to_from = ybus.branch_admittance_to_from
    end

    return Ybus(
        data,
        adjacency_data,
        (bus_ax, bus_ax),
        (bus_lookup, bus_lookup),
        subnetwork_axes,
        arc_subnetwork_axis,
        nr,
        branch_admittance_from_to,
        branch_admittance_to_from,
    )
end

function _make_subnetwork_axes(ybus, bus_numbers_to_remove, arcs_to_remove)
    subnetwork_axes = deepcopy(ybus.subnetwork_axes)
    arc_subnetwork_axis = deepcopy(ybus.arc_subnetwork_axis)
    for (k, values) in subnetwork_axes
        if k in bus_numbers_to_remove
            @warn "Reference bus removed during reduction; assigning arbitrary reference bus."
            axis_1, axis_2 = pop!(subnetwork_axes, k)
            new_ref_bus = pop!(axis_1)
            pop!(axis_2)
            subnetwork_axes[new_ref_bus] = (axis_1, axis_2)
            # If a reference bus key is reduced, change the arc subnetwork axis key as well: 
            arc_subnetwork_axis[new_ref_bus] = pop!(arc_subnetwork_axis, k)
        end
    end
    for (k, values) in subnetwork_axes
        new_values = setdiff(values[1], bus_numbers_to_remove)
        subnetwork_axes[k] = (new_values, new_values)
    end
    for (k, values) in arc_subnetwork_axis
        arc_subnetwork_axis[k] = setdiff(values, arcs_to_remove)
    end
    return subnetwork_axes, arc_subnetwork_axis
end

function _add_series_branches_to_ybus!(
    data::SparseArrays.SparseMatrixCSC{ComplexF32, Int64},
    bus_lookup::Dict{Int, Int},
    yft::Nothing,
    ytf::Nothing,
    series_branch_map::Dict{Tuple{Int, Int}, Vector{Any}},
    nrd,
)
    for (equivalent_arc, series_map_entry) in series_branch_map
        ordered_bus_numbers, segment_orientations =
            _get_chain_data(equivalent_arc, series_map_entry, nrd)
        ordered_bus_indices = [bus_lookup[x] for x in ordered_bus_numbers]
        equivalent_arc_indices = (ordered_bus_indices[1], ordered_bus_indices[end])
        ybus_isolated_d2_chain = _build_chain_ybus(series_map_entry, segment_orientations)
        ybus_boundary_isolated_d2_chain = _reduce_internal_nodes(ybus_isolated_d2_chain)
        _apply_d2_chain_ybus!(
            data,
            ybus_isolated_d2_chain,
            ybus_boundary_isolated_d2_chain,
            equivalent_arc_indices,
        )
    end
    return yft, ytf
end

function _add_series_branches_to_ybus!(
    data::SparseArrays.SparseMatrixCSC{ComplexF32, Int64},
    bus_lookup::Dict{Int, Int},
    yft::ArcAdmittanceMatrix,
    ytf::ArcAdmittanceMatrix,
    series_branch_map::Dict{Tuple{Int, Int}, Vector{Any}},
    nrd::NetworkReductionData,
)
    arc_lookup = get_arc_lookup(yft)
    arc_axis = get_arc_axis(yft)
    I_yft, J_yft, V_yft = SparseArrays.findnz(yft.data)
    I_ytf, J_ytf, V_ytf = SparseArrays.findnz(ytf.data)
    row_ix = size(yft)[1] + 1
    n_buses = size(yft)[2]
    for (equivalent_arc, series_map_entry) in series_branch_map
        ordered_bus_numbers, segment_orientations =
            _get_chain_data(equivalent_arc, series_map_entry, nrd)
        ordered_bus_indices = [bus_lookup[x] for x in ordered_bus_numbers]
        equivalent_arc_indices = (ordered_bus_indices[1], ordered_bus_indices[end])
        ybus_isolated_d2_chain = _build_chain_ybus(series_map_entry, segment_orientations)
        ybus_boundary_isolated_d2_chain = _reduce_internal_nodes(ybus_isolated_d2_chain)
        _apply_d2_chain_ybus!(
            data,
            ybus_isolated_d2_chain,
            ybus_boundary_isolated_d2_chain,
            equivalent_arc_indices,
        )

        push!(arc_axis, equivalent_arc)
        push!(I_yft, row_ix)
        push!(I_yft, row_ix)
        push!(J_yft, equivalent_arc_indices[1])
        push!(J_yft, equivalent_arc_indices[2])
        push!(V_yft, ybus_boundary_isolated_d2_chain[1, 1])
        push!(V_yft, ybus_boundary_isolated_d2_chain[1, 2])
        push!(I_ytf, row_ix)
        push!(I_ytf, row_ix)
        push!(J_ytf, equivalent_arc_indices[2])
        push!(J_ytf, equivalent_arc_indices[1])
        push!(V_ytf, ybus_boundary_isolated_d2_chain[2, 2])
        push!(V_ytf, ybus_boundary_isolated_d2_chain[2, 1])
        row_ix += 1
    end
    yft_data = SparseArrays.sparse(I_yft, J_yft, V_yft, row_ix - 1, n_buses)
    ytf_data = SparseArrays.sparse(I_ytf, J_ytf, V_ytf, row_ix - 1, n_buses)

    branch_admittance_from_to = ArcAdmittanceMatrix(
        yft_data,
        (arc_axis, get_bus_axis(yft)),
        (arc_lookup, get_bus_lookup(yft)),
        nrd,
        :FromTo,
    )
    branch_admittance_to_from = ArcAdmittanceMatrix(
        ytf_data,
        (arc_axis, get_bus_axis(ytf)),
        (arc_lookup, get_bus_lookup(ytf)),
        nrd,
        :ToFrom,
    )
    return branch_admittance_from_to, branch_admittance_to_from
end

function _apply_d2_chain_ybus!(
    ybus_full::SparseArrays.SparseMatrixCSC{ComplexF32, Int64},
    ybus_chain::Matrix{ComplexF32},
    ybus_chain_reduced::Matrix{ComplexF32},
    equivalent_arc_indices::Tuple{Int, Int},
)
    equivalent_from_index, equivalent_to_index = equivalent_arc_indices
    from_bus_entry_difference = ybus_chain_reduced[1, 1] - ybus_chain[1, 1]
    from_to_bus_entry = ybus_chain_reduced[1, 2]
    to_from_bus_entry = ybus_chain_reduced[2, 1]
    to_bus_entry_difference = ybus_chain_reduced[end, end] - ybus_chain[end, end]

    ybus_full[equivalent_from_index, equivalent_from_index] += from_bus_entry_difference
    ybus_full[equivalent_from_index, equivalent_to_index] = from_to_bus_entry
    ybus_full[equivalent_to_index, equivalent_from_index] = to_from_bus_entry
    ybus_full[equivalent_to_index, equivalent_to_index] += to_bus_entry_difference
    return
end

function _build_chain_ybus(series_chain::Vector{Any}, segment_orientations::Vector{Symbol})
    fb = Vector{Int}()
    tb = Vector{Int}()
    y11 = Vector{ComplexF32}()
    y12 = Vector{ComplexF32}()
    y21 = Vector{ComplexF32}()
    y22 = Vector{ComplexF32}()
    for (ix, segment) in enumerate(series_chain)
        add_segment_to_ybus!(
            segment,
            y11,
            y12,
            y21,
            y22,
            fb,
            tb,
            ix,
            segment_orientations[ix],
        )
    end
    return Matrix(
        SparseArrays.sparse(
            [fb; fb; tb; tb],  # row indices
            [fb; tb; fb; tb],  # column indices
            [y11; y12; y21; y22],  # values
        ),
    )
end

function add_segment_to_ybus!(
    segment::Union{PSY.Branch, Tuple{PSY.ThreeWindingTransformer, Int}},
    y11::Vector{ComplexF32},
    y12::Vector{ComplexF32},
    y21::Vector{ComplexF32},
    y22::Vector{ComplexF32},
    fb::Vector{Int},
    tb::Vector{Int},
    ix::Int,
    segment_orientation::Symbol,
)
    (Y11, Y12, Y21, Y22) = ybus_branch_entries(segment)
    push!(fb, ix)
    push!(tb, ix + 1)
    if segment_orientation == :FromTo
        push!(y11, Y11)
        push!(y12, Y12)
        push!(y21, Y21)
        push!(y22, Y22)
    elseif segment_orientation == :ToFrom
        push!(y11, Y22)
        push!(y12, Y21)
        push!(y21, Y12)
        push!(y22, Y11)
    else
        error("Invalid segment orientation $(segment_orientation)")
    end
end
function add_segment_to_ybus!(
    segment::Set{PSY.Branch},
    y11::Vector{ComplexF32},
    y12::Vector{ComplexF32},
    y21::Vector{ComplexF32},
    y22::Vector{ComplexF32},
    fb::Vector{Int},
    tb::Vector{Int},
    ix::Int,
    segment_orientation::Symbol,
)
    for branch in segment
        add_segment_to_ybus!(branch, y11, y12, y21, y22, fb, tb, ix, segment_orientation)
    end
end

function _get_chain_data(
    equivalent_arc::Tuple{Int, Int},
    series_map_entry::Vector{Any},
    nrd::NetworkReductionData,
)
    ordered_bus_numbers = [equivalent_arc[1]]
    segment_orientation = Vector{Symbol}()
    for segment in series_map_entry
        arc_tuple = get_arc_tuple(segment, nrd)
        if arc_tuple[1] in ordered_bus_numbers
            push!(ordered_bus_numbers, arc_tuple[2])
            push!(segment_orientation, :FromTo)
        elseif arc_tuple[2] in ordered_bus_numbers
            push!(ordered_bus_numbers, arc_tuple[1])
            push!(segment_orientation, :ToFrom)
        else
            error("Found disconnected series chain")
        end
    end
    @assert ordered_bus_numbers[end] == equivalent_arc[2]
    return ordered_bus_numbers, segment_orientation
end

function _reduce_internal_nodes(Y::Matrix{ComplexF32})
    dim_Y = size(Y)[1]
    keep_ix = [1, dim_Y]
    eliminate_ix = collect(2:(dim_Y - 1))
    Y_kk = Y[keep_ix, keep_ix]
    Y_ee = Y[eliminate_ix, eliminate_ix]
    Y_ke = Y[keep_ix, eliminate_ix]
    Y_ek = Y[eliminate_ix, keep_ix]
    Y_reduced = Y_kk - Y_ke * (Y_ee \ Y_ek)
    return Y_reduced
end

function _remake_reverse_direct_branch_map!(nr::NetworkReductionData)
    reverse_direct_branch_map = Dict{PSY.Branch, Tuple{Int, Int}}()
    for (k, v) in nr.direct_branch_map
        reverse_direct_branch_map[v] = k
    end
    nr.reverse_direct_branch_map = reverse_direct_branch_map
    return
end
function _remake_reverse_parallel_branch_map!(nr::NetworkReductionData)
    reverse_parallel_branch_map = Dict{PSY.Branch, Tuple{Int, Int}}()
    for (k, v) in nr.parallel_branch_map
        for x in v
            reverse_parallel_branch_map[x] = k
        end
    end
    nr.reverse_parallel_branch_map = reverse_parallel_branch_map
    return
end
function _remake_reverse_series_branch_map!(nr::NetworkReductionData)
    reverse_series_branch_map = Dict{PSY.Branch, Tuple{Int, Int}}()
    for (k, v) in nr.series_branch_map
        for x in v
            reverse_series_branch_map[x] = k
        end
    end
    nr.reverse_series_branch_map = reverse_series_branch_map
    return
end

function _remake_reverse_transformer3W_map!(nr::NetworkReductionData)
    reverse_transformer3W_map =
        Dict{Tuple{PSY.ThreeWindingTransformer, Int}, Tuple{Int, Int}}()
    for (k, v) in nr.transformer3W_map
        reverse_transformer3W_map[v] = k
    end
    nr.reverse_transformer3W_map = reverse_transformer3W_map
    return
end

"""
Validates connectivity by checking that the number of subnetworks is 1 (fully connected network).
"""
function validate_connectivity(M::Ybus)
    sub_nets = find_subnetworks(M)
    return length(sub_nets) == 1
end

"""
Evaluates subnetworks by looking for the subsets of buses connected each other,
but not connected with buses of other subsets.
"""
function find_subnetworks(M::Ybus)
    bus_numbers = M.axes[2]
    return find_subnetworks(M.adjacency_data, bus_numbers)
end

function get_reduction(
    ybus::Ybus,
    sys::PSY.System,
    reduction::DegreeTwoReduction,
)
    A = AdjacencyMatrix(ybus)
    return get_reduction(A, sys, reduction)
end

function _validate_study_buses(
    ybus::Ybus,
    study_buses::Vector{Int},
)
    #TODO - improve building the vector/set of valid bus numbers
    valid_bus_numbers = Set{Int}()
    for (k, v) in get_network_reduction_data(ybus).bus_reduction_map
        push!(valid_bus_numbers, k)
        union!(valid_bus_numbers, v)
    end
    for b in study_buses
        b ∉ valid_bus_numbers &&
            throw(IS.DataFormatError("Study bus $b not found in system"))
    end
    slack_bus_numbers = get_ref_bus(ybus)
    subnetwork_axes = ybus.subnetwork_axes
    for axes in values(subnetwork_axes)
        subnetwork_bus_ax = axes[1]
        all_in = all(x -> x in Set(subnetwork_bus_ax), study_buses)
        none_in = all(x -> !(x in Set(subnetwork_bus_ax)), study_buses)
        if Set(subnetwork_bus_ax) == Set(study_buses)
            @warn "The study buses comprise an entire island; ward reduction will not modify this island and other islands will be eliminated"
        end
        if !(all_in || none_in)
            throw(
                IS.DataFormatError(
                    "All study_buses must occur in a single synchronously connected system.",
                ),
            )
        end
        for sb in slack_bus_numbers
            if sb in subnetwork_bus_ax && sb ∉ study_buses && !(none_in)
                throw(
                    IS.DataFormatError(
                        "Slack bus $sb must be included in the study buses for an area that is partially reduced",
                    ),
                )
            end
        end
    end

    return
end

function get_reduction(
    ybus::Ybus,
    ::PSY.System,
    reduction::WardReduction,
)
    study_buses = get_study_buses(reduction)
    _validate_study_buses(ybus, study_buses)
    bus_lookup = get_bus_lookup(ybus)
    bus_axis = get_bus_axis(ybus)
    A = IncidenceMatrix(ybus)
    boundary_buses = Set{Int}()
    removed_arcs = Set{Tuple{Int, Int}}()
    for arc in get_arc_axis(A)
        #Deterimine boundary buses:
        if (arc[1] ∈ study_buses) && (arc[2] ∉ study_buses)
            push!(boundary_buses, arc[1])
        elseif (arc[1] ∉ study_buses) && (arc[2] ∈ study_buses)
            push!(boundary_buses, arc[2])
        end
        #Determine arcs outside of study area
        if !(arc[1] ∈ study_buses && arc[2] ∈ study_buses)
            push!(removed_arcs, arc)
        end
    end
    bus_reduction_map, reverse_bus_search_map, added_branch_map, added_admittance_map =
        get_ward_reduction(
            ybus.data,
            bus_lookup,
            bus_axis,
            boundary_buses,
            Set(get_ref_bus(ybus)),
            study_buses,
        )

    return NetworkReductionData(;
        bus_reduction_map = bus_reduction_map,
        reverse_bus_search_map = reverse_bus_search_map,
        removed_arcs = removed_arcs,
        added_branch_map = added_branch_map,
        added_admittance_map = added_admittance_map,
        reductions = ReductionContainer(; ward_reduction = reduction),
    )
end
