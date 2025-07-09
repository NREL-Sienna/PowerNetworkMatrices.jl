"""
Nodal admittance matrix (Ybus) is an N x N matrix describing a power system with N buses. It represents the nodal admittance of the buses in a power system.

The Ybus Struct is indexed using the Bus Numbers, no need for them to be sequential

The fields yft and ytf are the branch admittance matrices for the from-to and to-from branch admittances respectively. The rows correspond to branches and the columns to buses.
The matrix columns are mapped to buses using fb, tb arrays of the matrix columns that correspond to the `from` and `to` buses.
Using yft, ytf, and the voltage vector, the branch currents and power flows can be calculated.
"""
struct Ybus{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{ComplexF64}
    data::SparseArrays.SparseMatrixCSC{ComplexF64, Int}
    adjacency_data::SparseArrays.SparseMatrixCSC{Int8, Int}
    ref_bus_numbers::Set{Int}
    axes::Ax
    lookup::L
    subnetworks::Dict{Int, Set{Int}}
    network_reduction::NetworkReduction
    yft::Union{SparseArrays.SparseMatrixCSC{ComplexF64, Int}, Nothing}
    ytf::Union{SparseArrays.SparseMatrixCSC{ComplexF64, Int}, Nothing}
    fb::Union{Vector{Int64}, Nothing}
    tb::Union{Vector{Int64}, Nothing}
end

get_network_reduction(y::Ybus) = y.network_reduction

function add_to_branch_maps!(nr::NetworkReduction, arc::PSY.Arc, br::PSY.Branch)
    direct_branch_map = get_direct_branch_map(nr)
    reverse_direct_branch_map = get_reverse_direct_branch_map(nr)
    parallel_branch_map = get_parallel_branch_map(nr)
    reverse_parallel_branch_map = get_reverse_parallel_branch_map(nr)
    arc_tuple = get_arc_tuple(arc)
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
    nr::NetworkReduction,
    primary_star_arc::PSY.Arc,
    secondary_star_arc::PSY.Arc,
    tertiary_star_arc::PSY.Arc,
    br::Union{PSY.Transformer3W, PSY.PhaseShiftingTransformer3W},
)
    transformer3W_map = get_transformer3W_map(nr)
    reverse_transformer3W_map = get_reverse_transformer3W_map(nr)
    if PSY.get_available_primary(br)
        primary_star_arc_tuple = get_arc_tuple(primary_star_arc)
        transformer3W_map[primary_star_arc_tuple] = (br, 1)
        reverse_transformer3W_map[(br, 1)] = primary_star_arc_tuple
    end
    if PSY.get_available_secondary(br)
        secondary_star_arc_tuple = get_arc_tuple(secondary_star_arc)
        transformer3W_map[secondary_star_arc_tuple] = (br, 2)
        reverse_transformer3W_map[(br, 2)] = secondary_star_arc_tuple
    end
    if PSY.get_available_tertiary(br)
        tertiary_star_arc_tuple = get_arc_tuple(tertiary_star_arc)
        transformer3W_map[tertiary_star_arc_tuple] = (br, 3)
        reverse_transformer3W_map[(br, 3)] = tertiary_star_arc_tuple
    end
    return
end

function _ybus!(
    y11::Vector{ComplexF64},
    y12::Vector{ComplexF64},
    y21::Vector{ComplexF64},
    y22::Vector{ComplexF64},
    br::PSY.ACTransmission,
    num_bus::Dict{Int, Int},
    branch_ix::Int64,
    fb::Vector{Int64},
    tb::Vector{Int64},
    nr::NetworkReduction,
    adj::SparseArrays.SparseMatrixCSC{Int8, Int},
)
    arc = PSY.get_arc(br)
    add_to_branch_maps!(nr, arc, br)
    bus_from_no, bus_to_no = get_bus_indices(arc, num_bus, nr)
    adj[bus_from_no, bus_to_no] = 1
    adj[bus_to_no, bus_from_no] = -1
    fb[branch_ix] = bus_from_no
    tb[branch_ix] = bus_to_no
    Y_l = (1 / (PSY.get_r(br) + PSY.get_x(br) * 1im))
    Y11 = Y_l + (1im * PSY.get_b(br).from)
    if !isfinite(Y11) || !isfinite(Y_l)
        error(
            "Data in $(PSY.get_name(br)) is incorrect. r = $(PSY.get_r(br)), x = $(PSY.get_x(br))",
        )
    end
    y11[branch_ix] = Y11
    Y12 = -Y_l
    y12[branch_ix] = Y12
    Y21 = Y12
    y21[branch_ix] = Y21
    Y22 = Y_l + (1im * PSY.get_b(br).to)
    y22[branch_ix] = Y22
    return
end

function _ybus!(
    y11::Vector{ComplexF64},
    y12::Vector{ComplexF64},
    y21::Vector{ComplexF64},
    y22::Vector{ComplexF64},
    br::PSY.DiscreteControlledACBranch,
    num_bus::Dict{Int, Int},
    branch_ix::Int64,
    fb::Vector{Int64},
    tb::Vector{Int64},
    nr::NetworkReduction,
    adj::SparseArrays.SparseMatrixCSC{Int8, Int},
)
    arc = PSY.get_arc(br)
    add_to_branch_maps!(nr, arc, br)
    bus_from_no, bus_to_no = get_bus_indices(arc, num_bus, nr)
    adj[bus_from_no, bus_to_no] = 1
    adj[bus_to_no, bus_from_no] = -1
    fb[branch_ix] = bus_from_no
    tb[branch_ix] = bus_to_no
    Y_l = (1 / (PSY.get_r(br) + PSY.get_x(br) * 1im))
    Y11 = Y_l
    if !isfinite(Y11) || !isfinite(Y_l)
        error(
            "Data in $(PSY.get_name(br)) is incorrect. r = $(PSY.get_r(br)), x = $(PSY.get_x(br))",
        )
    end
    y11[branch_ix] = Y11
    Y12 = -Y_l
    y12[branch_ix] = Y12
    Y21 = Y12
    y21[branch_ix] = Y21
    Y22 = Y_l
    y22[branch_ix] = Y22
    return
end

function _ybus!(
    y11::Vector{ComplexF64},
    y12::Vector{ComplexF64},
    y21::Vector{ComplexF64},
    y22::Vector{ComplexF64},
    br::PSY.DynamicBranch,
    num_bus::Dict{Int, Int},
    reverse_bus_search_map::Dict{Int, Int},
    branch_ix::Int64,
    fb::Vector{Int64},
    tb::Vector{Int64},
    nr::NetworkReduction,
    adj::SparseArrays.SparseMatrixCSC{Int8, Int},
)
    _ybus!(
        y11,
        y12,
        y21,
        y22,
        br.branch,
        num_bus,
        reverse_bus_search_map,
        branch_ix,
        fb,
        tb,
        nr,
        adj,
    )
    return
end

function _ybus!(
    y11::Vector{ComplexF64},
    y12::Vector{ComplexF64},
    y21::Vector{ComplexF64},
    y22::Vector{ComplexF64},
    br::PSY.Transformer2W,
    num_bus::Dict{Int, Int},
    branch_ix::Int64,
    fb::Vector{Int64},
    tb::Vector{Int64},
    nr::NetworkReduction,
    adj::SparseArrays.SparseMatrixCSC{Int8, Int},
)
    arc = PSY.get_arc(br)
    add_to_branch_maps!(nr, arc, br)
    bus_from_no, bus_to_no = get_bus_indices(arc, num_bus, nr)
    adj[bus_from_no, bus_to_no] = 1
    adj[bus_to_no, bus_from_no] = -1
    fb[branch_ix] = bus_from_no
    tb[branch_ix] = bus_to_no
    Y_t = 1 / (PSY.get_r(br) + PSY.get_x(br) * 1im)
    Y11 = Y_t
    y_shunt = PSY.get_primary_shunt(br)
    if !isfinite(Y11) || !isfinite(Y_t) || !isfinite(y_shunt)
        error(
            "Data in $(PSY.get_name(br)) is incorrect. r = $(PSY.get_r(br)), x = $(PSY.get_x(br))",
        )
    end

    y11[branch_ix] = Y11 + y_shunt
    Y12 = -Y_t
    y12[branch_ix] = Y12
    Y21 = Y12
    y21[branch_ix] = Y21
    Y22 = Y_t
    y22[branch_ix] = Y22
    return
end

function _ybus!(
    y11::Vector{ComplexF64},
    y12::Vector{ComplexF64},
    y21::Vector{ComplexF64},
    y22::Vector{ComplexF64},
    br::PSY.Transformer3W,
    num_bus::Dict{Int, Int},
    offset_ix::Int64,
    fb::Vector{Int64},
    tb::Vector{Int64},
    ix::Int64,
    nr::NetworkReduction,
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
    y_shunt = PSY.get_g(br) + im * PSY.get_b(br)
    if primary_available
        primary_ix, star_ix = get_bus_indices(primary_star_arc, num_bus, nr)
        adj[primary_ix, star_ix] = 1
        adj[star_ix, primary_ix] = -1
        fb[offset_ix + ix + n_entries] = primary_ix
        tb[offset_ix + ix + n_entries] = star_ix
        c = 1 / PSY.get_primary_turns_ratio(br)
        Y_t = 1 / (PSY.get_r_primary(br) + PSY.get_x_primary(br) * 1im)
        Y11 = (Y_t * c^2)
        if !isfinite(Y11) || !isfinite(y_shunt)
            error(
                "Data in $(PSY.get_name(br)) is incorrect.
                r_p = $(PSY.get_r_primary(br)), x_p = $(PSY.get_x_primary(br))",
            )
        end
        y11[offset_ix + ix + n_entries] = Y11 + y_shunt
        Y12 = (-Y_t * c)
        y12[offset_ix + ix + n_entries] = Y12
        Y21 = Y12
        y21[offset_ix + ix + n_entries] = Y21
        Y22 = Y_t
        y22[offset_ix + ix + n_entries] = Y22
        n_entries += 1
    end
    if secondary_available
        secondary_ix, star_ix = get_bus_indices(secondary_star_arc, num_bus, nr)
        adj[secondary_ix, star_ix] = 1
        adj[star_ix, secondary_ix] = -1
        fb[offset_ix + ix + n_entries] = secondary_ix
        tb[offset_ix + ix + n_entries] = star_ix
        c = 1 / PSY.get_secondary_turns_ratio(br)
        Y_t = 1 / (PSY.get_r_secondary(br) + PSY.get_x_secondary(br) * 1im)
        Y11 = (Y_t * c^2)
        if !isfinite(Y11) || !isfinite(y_shunt)
            error(
                "Data in $(PSY.get_name(br)) is incorrect.
                r_p = $(PSY.get_r_primary(br)), x_p = $(PSY.get_x_primary(br))",
            )
        end
        y11[offset_ix + ix + n_entries] = Y11
        Y12 = (-Y_t * c)
        y12[offset_ix + ix + n_entries] = Y12
        Y21 = Y12
        y21[offset_ix + ix + n_entries] = Y21
        Y22 = Y_t
        y22[offset_ix + ix + n_entries] = Y22
        n_entries += 1
    end
    if tertiary_available
        tertiary_ix, star_ix = get_bus_indices(tertiary_star_arc, num_bus, nr)
        adj[tertiary_ix, star_ix] = 1
        adj[star_ix, tertiary_ix] = -1
        fb[offset_ix + ix + n_entries] = tertiary_ix
        tb[offset_ix + ix + n_entries] = star_ix
        c = 1 / PSY.get_tertiary_turns_ratio(br)
        Y_t = 1 / (PSY.get_r_tertiary(br) + PSY.get_x_tertiary(br) * 1im)
        Y11 = (Y_t * c^2)
        if !isfinite(Y11) || !isfinite(y_shunt)
            error(
                "Data in $(PSY.get_name(br)) is incorrect.
                r_p = $(PSY.get_r_primary(br)), x_p = $(PSY.get_x_primary(br))",
            )
        end
        y11[offset_ix + ix + n_entries] = Y11
        Y12 = (-Y_t * c)
        y12[offset_ix + ix + n_entries] = Y12
        Y21 = Y12
        y21[offset_ix + ix + n_entries] = Y21
        Y22 = Y_t
        y22[offset_ix + ix + n_entries] = Y22
        n_entries += 1
    end
    return n_entries
end

function _ybus!(
    y11::Vector{ComplexF64},
    y12::Vector{ComplexF64},
    y21::Vector{ComplexF64},
    y22::Vector{ComplexF64},
    br::PSY.PhaseShiftingTransformer3W,
    num_bus::Dict{Int, Int},
    offset_ix::Int64,
    fb::Vector{Int64},
    tb::Vector{Int64},
    ix::Int64,
    nr::NetworkReduction,
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
    y_shunt = PSY.get_g(br) + im * PSY.get_b(br)
    if primary_available
        primary_ix, star_ix = get_bus_indices(primary_star_arc, num_bus, nr)
        adj[primary_ix, star_ix] = 1
        adj[star_ix, primary_ix] = -1
        fb[offset_ix + ix + n_entries] = primary_ix
        tb[offset_ix + ix + n_entries] = star_ix
        Y_t = 1 / (PSY.get_r_primary(br) + PSY.get_x_primary(br) * 1im)
        tap = (PSY.get_primary_turns_ratio(br) * exp(PSY.get_α_primary(br) * 1im))
        c_tap = (PSY.get_primary_turns_ratio(br) * exp(-1 * PSY.get_α_primary(br) * 1im))
        Y11 = (Y_t / abs(tap)^2)
        if !isfinite(Y11) || !isfinite(y_shunt)
            error(
                "Data in $(PSY.get_name(br)) is incorrect.
                r_p = $(PSY.get_r_primary(br)), x_p = $(PSY.get_x_primary(br))",
            )
        end
        y11[offset_ix + ix + n_entries] = Y11 + y_shunt
        Y12 = (-Y_t / c_tap)
        y12[offset_ix + ix + n_entries] = Y12
        Y21 = (-Y_t / tap)
        y21[offset_ix + ix + n_entries] = Y21
        Y22 = Y_t
        y22[offset_ix + ix + n_entries] = Y22
        n_entries += 1
    end
    if secondary_available
        secondary_ix, star_ix = get_bus_indices(secondary_star_arc, num_bus, nr)
        adj[secondary_ix, star_ix] = 1
        adj[star_ix, secondary_ix] = -1
        fb[offset_ix + ix + n_entries] = secondary_ix
        tb[offset_ix + ix + n_entries] = star_ix
        Y_t = 1 / (PSY.get_r_secondary(br) + PSY.get_x_secondary(br) * 1im)
        tap = (PSY.get_secondary_turns_ratio(br) * exp(PSY.get_α_secondary(br) * 1im))
        c_tap =
            (PSY.get_secondary_turns_ratio(br) * exp(-1 * PSY.get_α_secondary(br) * 1im))
        Y11 = (Y_t / abs(tap)^2)
        if !isfinite(Y11) || !isfinite(y_shunt)
            error(
                "Data in $(PSY.get_name(br)) is incorrect.
                r_p = $(PSY.get_r_primary(br)), x_p = $(PSY.get_x_primary(br))",
            )
        end
        y11[offset_ix + ix + n_entries] = Y11
        Y12 = (-Y_t / c_tap)
        y12[offset_ix + ix + n_entries] = Y12
        Y21 = (-Y_t / tap)
        y21[offset_ix + ix + n_entries] = Y21
        Y22 = Y_t
        y22[offset_ix + ix + n_entries] = Y22
        n_entries += 1
    end
    if tertiary_available
        tertiary_ix, star_ix = get_bus_indices(tertiary_star_arc, num_bus, nr)
        adj[tertiary_ix, star_ix] = 1
        adj[star_ix, tertiary_ix] = -1
        fb[offset_ix + ix + n_entries] = tertiary_ix
        tb[offset_ix + ix + n_entries] = star_ix
        Y_t = 1 / (PSY.get_r_tertiary(br) + PSY.get_x_tertiary(br) * 1im)
        tap = (PSY.get_tertiary_turns_ratio(br) * exp(PSY.get_α_tertiary(br) * 1im))
        c_tap = (PSY.get_tertiary_turns_ratio(br) * exp(-1 * PSY.get_α_tertiary(br) * 1im))
        Y11 = (Y_t / abs(tap)^2)
        if !isfinite(Y11) || !isfinite(y_shunt)
            error(
                "Data in $(PSY.get_name(br)) is incorrect.
                r_p = $(PSY.get_r_primary(br)), x_p = $(PSY.get_x_primary(br))",
            )
        end
        y11[offset_ix + ix + n_entries] = Y11
        Y12 = (-Y_t / c_tap)
        y12[offset_ix + ix + n_entries] = Y12
        Y21 = (-Y_t / tap)
        y21[offset_ix + ix + n_entries] = Y21
        Y22 = Y_t
        y22[offset_ix + ix + n_entries] = Y22
        n_entries += 1
    end
    return n_entries
end

function _ybus!(
    y11::Vector{ComplexF64},
    y12::Vector{ComplexF64},
    y21::Vector{ComplexF64},
    y22::Vector{ComplexF64},
    br::PSY.TapTransformer,
    num_bus::Dict{Int, Int},
    branch_ix::Int64,
    fb::Vector{Int64},
    tb::Vector{Int64},
    nr::NetworkReduction,
    adj::SparseArrays.SparseMatrixCSC{Int8, Int},
)
    arc = PSY.get_arc(br)
    add_to_branch_maps!(nr, arc, br)
    bus_from_no, bus_to_no = get_bus_indices(arc, num_bus, nr)
    adj[bus_from_no, bus_to_no] = 1
    adj[bus_to_no, bus_from_no] = -1
    fb[branch_ix] = bus_from_no
    tb[branch_ix] = bus_to_no

    Y_t = 1 / (PSY.get_r(br) + PSY.get_x(br) * 1im)
    c = 1 / PSY.get_tap(br)
    y_shunt = PSY.get_primary_shunt(br)

    Y11 = (Y_t * c^2)
    y11[branch_ix] = Y11 + y_shunt
    Y12 = (-Y_t * c)
    if !isfinite(Y11) || !isfinite(Y12) || !isfinite(y_shunt)
        error(
            "Data in $(PSY.get_name(br)) is incorrect. r = $(PSY.get_r(br)), x = $(PSY.get_x(br))",
        )
    end

    y12[branch_ix] = Y12
    Y21 = Y12
    y21[branch_ix] = Y21
    Y22 = Y_t
    y22[branch_ix] = Y22
    return
end

function _ybus!(
    y11::Vector{ComplexF64},
    y12::Vector{ComplexF64},
    y21::Vector{ComplexF64},
    y22::Vector{ComplexF64},
    br::PSY.PhaseShiftingTransformer,
    num_bus::Dict{Int, Int},
    branch_ix::Int64,
    fb::Vector{Int64},
    tb::Vector{Int64},
    nr::NetworkReduction,
    adj::SparseArrays.SparseMatrixCSC{Int8, Int},
)
    arc = PSY.get_arc(br)
    add_to_branch_maps!(nr, arc, br)
    bus_from_no, bus_to_no = get_bus_indices(arc, num_bus, nr)
    adj[bus_from_no, bus_to_no] = 1
    adj[bus_to_no, bus_from_no] = -1
    fb[branch_ix] = bus_from_no
    tb[branch_ix] = bus_to_no
    Y_t = 1 / (PSY.get_r(br) + PSY.get_x(br) * 1im)
    tap = (PSY.get_tap(br) * exp(PSY.get_α(br) * 1im))
    c_tap = (PSY.get_tap(br) * exp(-1 * PSY.get_α(br) * 1im))
    y_shunt = PSY.get_primary_shunt(br)
    Y11 = (Y_t / abs(tap)^2)
    if !isfinite(Y11) || !isfinite(Y_t) || !isfinite(y_shunt * c_tap)
        error(
            "Data in $(PSY.get_name(br)) is incorrect. r = $(PSY.get_r(br)), x = $(PSY.get_x(br))",
        )
    end

    y11[branch_ix] = Y11 + y_shunt
    Y12 = (-Y_t / c_tap)
    y12[branch_ix] = Y12
    Y21 = (-Y_t / tap)
    y21[branch_ix] = Y21
    Y22 = Y_t
    y22[branch_ix] = Y22
    return
end

function _ybus!(
    ysh::Vector{ComplexF64},
    fa::PSY.FixedAdmittance,
    num_bus::Dict{Int, Int},
    fa_ix::Int64,
    sb::Vector{Int64},
    nr::NetworkReduction,
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
    ysh::Vector{ComplexF64},
    fa::PSY.SwitchedAdmittance,
    num_bus::Dict{Int, Int},
    fa_ix::Int64,
    sb::Vector{Int64},
    nr::NetworkReduction,
)
    bus_no = get_bus_index(fa, num_bus, nr)
    sb[fa_ix] = bus_no
    ysh[fa_ix] = 0.0
    return
end

function _ybus!(
    ysh::Vector{ComplexF64},
    fa::PSY.StandardLoad,
    num_bus::Dict{Int, Int},
    fa_ix::Int64,
    sb::Vector{Int64},
    nr::NetworkReduction,
)
    bus_no = get_bus_index(fa, num_bus, nr)
    Y = PSY.get_impedance_active_power(fa) + im * PSY.get_impedance_reactive_power(fa)
    sb[fa_ix] = bus_no
    if !isfinite(Y)
        error(
            "Data in $(PSY.get_name(fa)) is incorrect. Y = $(Y)",
        )
    end
    ysh[fa_ix] = Y
    return
end

function _ybus!(
    y11::Vector{ComplexF64},
    y12::Vector{ComplexF64},
    y21::Vector{ComplexF64},
    y22::Vector{ComplexF64},
    br::PSY.DiscreteControlledACBranch,
    num_bus::Dict{Int, Int},
    branch_ix::Int64,
    fb::Vector{Int64},
    tb::Vector{Int64},
    nr::NetworkReduction,
)
    arc = PSY.get_arc(br)
    bus_from_no, bus_to_no = get_bus_indices(arc, num_bus, nr)
    adj[bus_from_no, bus_to_no] = 1
    adj[bus_to_no, bus_from_no] = -1
    fb[branch_ix] = bus_from_no
    tb[branch_ix] = bus_to_no
    Y_l = (1 / (PSY.get_r(br) + PSY.get_x(br) * 1im))
    if !isfinite(Y_l)
        error(
            "Data in $(PSY.get_name(br)) is incorrect. r = $(PSY.get_r(br)), x = $(PSY.get_x(br))",
        )
    end
    y11[branch_ix] = Y_l
    y12[branch_ix] = -Y_l
    y21[branch_ix] = -Y_l
    y22[branch_ix] = Y_l
    return
end

function _buildybus!(
    network_reduction::NetworkReduction,
    adj::SparseArrays.SparseMatrixCSC{Int8, Int},
    branches,
    transformer_3w::Vector{Union{PSY.PhaseShiftingTransformer3W, PSY.Transformer3W}},
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
    fb = zeros(Int64, branchcount)
    tb = zeros(Int64, branchcount)
    sb = zeros(Int64, fa_count + sa_count + sl_count)

    y11 = zeros(ComplexF64, branchcount)
    y12 = zeros(ComplexF64, branchcount)
    y21 = zeros(ComplexF64, branchcount)
    y22 = zeros(ComplexF64, branchcount)
    ysh = zeros(ComplexF64, fa_count + sa_count + sl_count)

    for (ix, b) in enumerate(branches)
        if PSY.get_name(b) == "init"
            throw(DataFormatError("The data in Branch is invalid"))
        end
        _ybus!(y11, y12, y21, y22, b, num_bus, ix, fb, tb, network_reduction, adj)
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
            network_reduction,
            adj,
        )
        ix += n_entries
    end
    for (ix, fa) in enumerate([fixed_admittances; switched_admittances; standard_loads])
        _ybus!(ysh, fa, num_bus, ix, sb, network_reduction)
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

#= """
Builds a Ybus from a collection of buses and branches. The return is a Ybus Array indexed with the bus numbers and the branch names.

# Arguments
- `check_connectivity::Bool`: Checks connectivity of the network using Depth First Search (DFS)
"""
function Ybus(
    branches::Vector,
    buses::Vector{PSY.ACBus},
    transformer_3w::Vector{PSY.Transformer3W} = Vector{PSY.Transformer3W}(),
    fixed_admittances::Vector{PSY.FixedAdmittance} = Vector{PSY.FixedAdmittance}(),
    switched_admittances::Vector{PSY.SwitchedAdmittance} = Vector{PSY.SwitchedAdmittance}(),
    standard_loads::Vector{PSY.StandardLoad} = Vector{PSY.StandardLoad}();
    check_connectivity::Bool = true,
    make_branch_admittance_matrices::Bool = false,
    network_reduction = NetworkReduction(),
)
    bus_ax = PSY.get_number.(buses)
    axes = (bus_ax, bus_ax)
    bus_lookup = make_ax_ref(bus_ax)
    busnumber = length(buses)
    look_up = (bus_lookup, bus_lookup)
    adj = SparseArrays.spdiagm(ones(Int8, busnumber))

    y11, y12, y21, y22, ysh, fb, tb, sb =
        _buildybus!(
            network_reduction,
            adj,
            branches,
            transformer_3w,
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
    if check_connectivity && length(buses) > 1
        islands = find_subnetworks(ybus, bus_ax)
        length(islands) > 1 && throw(IS.DataFormatError("Network not connected"))
    end
    if make_branch_admittance_matrices
        yft = SparseArrays.sparse(
            [1:length(fb); 1:length(fb)],
            [fb; tb],
            [y11; y12],
            length(fb),
            length(buses),
        )
        ytf = SparseArrays.sparse(
            [1:length(tb); 1:length(tb)],
            [tb; fb],
            [y22; y21],
            length(tb),
            length(buses),
        )
    else
        yft = nothing
        ytf = nothing
        fb = nothing
        tb = nothing
    end
    return Ybus(ybus, adj, axes, look_up, network_reduction, yft, ytf, fb, tb)
end =#

"""
Builds a Ybus from the system. The return is a Ybus Array indexed with the bus numbers and the branch names.

# Arguments
- `check_connectivity::Bool`: Checks connectivity of the network
"""
function Ybus(
    sys::PSY.System;
    check_connectivity::Bool = true,
    make_branch_admittance_matrices::Bool = false,
    network_reductions::Vector{NetworkReductionTypes} = NetworkReductionTypes[],
    kwargs...,
)
    ref_bus_numbers = Set{Int}()
    nr = NetworkReduction()
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
                Union{PSY.Transformer3W, PSY.PhaseShiftingTransformer3W},
                sys,
            ),
        )
    fixed_admittances =
        collect(PSY.get_components(x -> PSY.get_available(x), PSY.FixedAdmittance, sys))
    switched_admittances =
        collect(PSY.get_components(x -> PSY.get_available(x), PSY.SwitchedAdmittance, sys))
    standard_loads =
        collect(PSY.get_components(x -> PSY.get_available(x), PSY.StandardLoad, sys))
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
        yft = SparseArrays.sparse(
            [1:length(fb); 1:length(fb)],
            [fb; tb],
            [y11; y12],
            length(fb),
            busnumber,
        )
        ytf = SparseArrays.sparse(
            [1:length(tb); 1:length(tb)],
            [tb; fb],
            [y22; y21],
            length(tb),
            busnumber,
        )
    else
        yft = nothing
        ytf = nothing
        fb = nothing
        tb = nothing
    end
    ybus = Ybus(
        ybus,
        adj,
        ref_bus_numbers,
        axes,
        lookup,
        Dict{Int, Set{Int}}(),
        nr,
        yft,
        ytf,
        fb,
        tb,
    )
    for nr in network_reductions
        ybus = build_reduced_ybus(ybus, sys, nr)
    end
    if check_connectivity && length(bus_lookup) > 1
        subnetworks = assign_reference_buses!(
            find_subnetworks(ybus.data, ybus.axes[1]),
            ybus.ref_bus_numbers,
        )
        if length(subnetworks) > 1
            @warn "More than one island found; Network is not connected"
        end
        return Ybus(
            ybus.data,
            ybus.adjacency_data,
            ybus.ref_bus_numbers,
            ybus.axes,
            ybus.lookup,
            subnetworks,
            ybus.network_reduction,
            ybus.yft,
            ybus.ytf,
            ybus.fb,
            ybus.tb,
        )
    else
        return ybus
    end
end

function _goderya(ybus::SparseArrays.SparseMatrixCSC)
    node_count = size(ybus)[1]
    max_I = node_count^2
    I, J, val = SparseArrays.findnz(ybus)
    T = SparseArrays.sparse(I, J, ones(Int, length(val)))
    T_ = T * T
    for n in 1:(node_count - 1)
        I, _, _ = SparseArrays.findnz(T_)
        if length(I) == max_I
            @info "The System has no islands"
            break
        elseif length(I) < max_I
            temp = T_ * T
            I_temp, _, _ = SparseArrays.findnz(temp)
            if all(I_temp == I)
                @warn "The system contains islands" maxlog = 1
            end
            T_ = temp
        else
            @assert false
        end
        #@assert n < node_count - 1
    end
    return I
end

function validate_connectivity(
    M,
    nodes::Vector{PSY.ACBus},
    bus_lookup::Dict{Int64, Int64};
    connectivity_method::Function = goderya_connectivity,
)
    connected = connectivity_method(M, nodes, bus_lookup)
    return connected
end

function goderya_connectivity(M, nodes::Vector{PSY.ACBus}, bus_lookup::Dict{Int64, Int64})
    @info "Validating connectivity with Goderya algorithm"
    length(nodes) > 15_000 &&
        @warn "The Goderya algorithm is memory intensive on large networks and may not scale well, try `connectivity_method = dfs_connectivity"

    I = _goderya(M)

    node_count = length(nodes)
    connections = Dict([i => count(x -> x == i, I) for i in Set(I)])

    if length(Set(I)) == node_count
        connected = true
        if any(values(connections) .!= node_count)
            cc = Set(values(connections))
            @warn "Network has at least $(length(cc)) connected components with $cc nodes"
            connected = false
        end
    else
        disconnected_nodes = PSY.get_name.(nodes[setdiff(values(bus_lookup), I)])
        @warn "Principal connected component does not contain:" disconnected_nodes
        connected = false
    end
    return connected
end

"""
Finds the set of bus numbers that belong to each connected component in the System
"""
# this function extends the PowerModels.jl implementation to accept a System
function find_connected_components(sys::PSY.System)
    a = Adjacency(sys; check_connectivity = false)
    return find_connected_components(a.data, a.lookup[1])
end

# this function extends the PowerModels.jl implementation to accept an adjacency matrix and bus lookup
function find_connected_components(M, bus_lookup::Dict{Int64, Int64})
    pm_buses = Dict([i => Dict("bus_type" => 1, "bus_i" => b) for (i, b) in bus_lookup])

    arcs = findall((LinearAlgebra.UpperTriangular(M) - LinearAlgebra.I) .!= 0)
    pm_branches = Dict([
        i => Dict("f_bus" => a[1], "t_bus" => a[2], "br_status" => 1) for
        (i, a) in enumerate(arcs)
    ],)

    data = Dict("bus" => pm_buses, "branch" => pm_branches)
    cc = PSY.calc_connected_components(data)
    bus_decode = Dict(value => key for (key, value) in bus_lookup)
    connected_components = Vector{Set{Int64}}()
    for c in cc
        push!(connected_components, Set([bus_decode[b] for b in c]))
    end
    return Set(connected_components)
end

function dfs_connectivity(M, ::Vector{PSY.ACBus}, bus_lookup::Dict{Int64, Int64})
    @info "Validating connectivity with depth first search (network traversal)"
    cc = find_connected_components(M, bus_lookup)
    if length(cc) != 1
        @warn "Network has at least $(length(cc)) connected components with $(length.(cc)) nodes"
        connected = false
    else
        connected = true
    end
    return connected
end

function build_reduced_ybus(ybus::Ybus, sys::PSY.System, reduction::NetworkReductionTypes)
    A = IncidenceMatrix(ybus)
    nr = get_reduction(A, sys, Val(reduction))
    return _apply_reduction(ybus, nr)
end

function _apply_reduction(ybus::Ybus, nr_new::NetworkReduction)
    remake_reverse_direct_branch_map = false
    remake_reverse_parallel_branch_map = false
    remake_reverse_series_branch_map = false
    remake_reverse_transformer3W_map = false
    data = get_data(ybus)
    adjacency_data = ybus.adjacency_data    #TODO add getters 
    lookup = get_lookup(ybus)
    bus_lookup = lookup[1]
    nr = ybus.network_reduction
    bus_numbers_to_remove = Vector{Int}()
    for (k, v) in nr_new.reverse_bus_search_map
        nr.reverse_bus_search_map[k] = v
        push!(bus_numbers_to_remove, k)
    end
    for x in nr_new.removed_buses
        push!(nr.removed_buses, x)
        push!(bus_numbers_to_remove, x)
    end
    for (k, v) in nr_new.bus_reduction_map
        if haskey(nr.bus_reduction_map, k)
            union!(nr.bus_reduction_map[k], nr_new.bus_reduction_map[k])
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
    # TODO - order matters; can we depend on ordered branches from the degree_two_reduction algorithm or need to compute here from arbitrary set of branches
    for (k, v) in nr_new.series_branch_map
        nr.series_branch_map[k] = v
        total_series_impedance_from_to = 0.0
        total_series_impedance_to_from = 0.0
        total_shunt_admittance = 0.0
        for x in v
            fr_bus_no = x.arc.from.number
            to_bus_no = x.arc.to.number
            total_series_impedance_from_to +=
                1 / data[bus_lookup[fr_bus_no], bus_lookup[to_bus_no]]
            total_series_impedance_to_from +=
                1 / data[bus_lookup[to_bus_no], bus_lookup[fr_bus_no]]
            if fr_bus_no ∉ k
                total_shunt_admittance += data[bus_lookup[fr_bus_no], bus_lookup[fr_bus_no]]
            end
            if to_bus_no ∉ k
                total_shunt_admittance += data[bus_lookup[to_bus_no], bus_lookup[to_bus_no]]
            end
        end
        data[bus_lookup[k[1]], bus_lookup[k[2]]] += 1 / total_series_impedance_from_to
        data[bus_lookup[k[2]], bus_lookup[k[1]]] += 1 / total_series_impedance_to_from
        data[bus_lookup[k[1]], bus_lookup[k[1]]] += total_shunt_admittance / 2
        data[bus_lookup[k[2]], bus_lookup[k[2]]] += total_shunt_admittance / 2
    end

    remake_reverse_direct_branch_map && _remake_reverse_direct_branch_map!(nr)
    remake_reverse_parallel_branch_map && _remake_reverse_parallel_branch_map!(nr)
    remake_reverse_series_branch_map && _remake_reverse_series_branch_map!(nr)
    remake_reverse_transformer3W_map && _remake_reverse_transformer3W_map!(nr)

    push!(nr.reduction_type, nr_new.reduction_type[1])
    bus_ax = setdiff(ybus.axes[1], bus_numbers_to_remove)
    bus_lookup = make_ax_ref(bus_ax)
    bus_ix = [ybus.lookup[1][x] for x in bus_ax]
    adjacency_data = adjacency_data[bus_ix, bus_ix]
    data = data[bus_ix, bus_ix]
    #TODO - doesn't yet account for reduction in yft, ytf, fb, tb
    return Ybus(
        data,
        adjacency_data,
        setdiff(ybus.ref_bus_numbers, bus_numbers_to_remove),
        (bus_ax, bus_ax),
        (bus_lookup, bus_lookup),
        Dict{Int, Set{Int}}(),
        nr,
        ybus.yft,
        ybus.ytf,
        ybus.fb,
        ybus.tb,
    )
end

function _remake_reverse_direct_branch_map!(nr::NetworkReduction)
    reverse_direct_branch_map = Dict{PSY.Branch, Tuple{Int, Int}}()
    for (k, v) in nr.direct_branch_map
        reverse_direct_branch_map[v] = k
    end
    nr.reverse_direct_branch_map = reverse_direct_branch_map
end
function _remake_reverse_parallel_branch_map!(nr::NetworkReduction)
    reverse_parallel_branch_map = Dict{PSY.Branch, Tuple{Int, Int}}()
    for (k, v) in nr.parallel_branch_map
        for x in v
            reverse_parallel_branch_map[x] = k
        end
    end
    nr.reverse_parallel_branch_map = reverse_parallel_branch_map
end
function _remake_reverse_series_branch_map!(nr::NetworkReduction)
    reverse_series_branch_map = Dict{PSY.Branch, Tuple{Int, Int}}()
    for (k, v) in nr.series_branch_map
        for x in v
            reverse_series_branch_map[x] = k
        end
    end
    nr.reverse_series_branch_map = reverse_series_branch_map
end
function _remake_reverse_transformer3W_map!(nr::NetworkReduction)
    reverse_transformer3W_map = Dict{Tuple{PSY.Transformer3W, Int}, Tuple{Int, Int}}()
    for (k, v) in nr.transformer3W_map
        reverse_transformer3W_map[v] = k
    end
    nr.reverse_transformer3W_map = reverse_transformer3W_map
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
