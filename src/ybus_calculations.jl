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
    ref_bus_numbers::Set{Int}
    axes::Ax
    lookup::L
    subnetworks::Dict{Int, Set{Int}}
    network_reduction_data::NetworkReductionData
    yft::Union{SparseArrays.SparseMatrixCSC{ComplexF32, Int}, Nothing}
    ytf::Union{SparseArrays.SparseMatrixCSC{ComplexF32, Int}, Nothing}
    fb::Union{Vector{Int}, Nothing}
    tb::Union{Vector{Int}, Nothing}
end

get_network_reduction_data(y::Ybus) = y.network_reduction_data

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
    y11::Vector{ComplexF32},
    y12::Vector{ComplexF32},
    y21::Vector{ComplexF32},
    y22::Vector{ComplexF32},
    br::PSY.DiscreteControlledACBranch,
    num_bus::Dict{Int, Int},
    branch_ix::Int,
    fb::Vector{Int},
    tb::Vector{Int},
    nr::NetworkReductionData,
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
    y11::Vector{ComplexF32},
    y12::Vector{ComplexF32},
    y21::Vector{ComplexF32},
    y22::Vector{ComplexF32},
    br::PSY.DynamicBranch,
    num_bus::Dict{Int, Int},
    reverse_bus_search_map::Dict{Int, Int},
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
    y11::Vector{ComplexF32},
    y12::Vector{ComplexF32},
    y21::Vector{ComplexF32},
    y22::Vector{ComplexF32},
    br::PSY.Transformer2W,
    num_bus::Dict{Int, Int},
    branch_ix::Int,
    fb::Vector{Int},
    tb::Vector{Int},
    nr::NetworkReductionData,
    adj::SparseArrays.SparseMatrixCSC{Int8, Int},
)
    arc = PSY.get_arc(br)
    add_to_branch_maps!(nr, arc, br)
    bus_from_no, bus_to_no = get_bus_indices(arc, num_bus, nr)
    Y_t = 1 / (PSY.get_r(br) + PSY.get_x(br) * 1im)
    Y11 = Y_t
    y_shunt = PSY.get_primary_shunt(br)
    if !isfinite(Y11) || !isfinite(Y_t) || !isfinite(y_shunt)
        error(
            "Data in $(PSY.get_name(br)) is incorrect. r = $(PSY.get_r(br)), x = $(PSY.get_x(br))",
        )
    end

    adj[bus_from_no, bus_to_no] = 1
    adj[bus_to_no, bus_from_no] = -1
    fb[branch_ix] = bus_from_no
    tb[branch_ix] = bus_to_no

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
    y11::Vector{ComplexF32},
    y12::Vector{ComplexF32},
    y21::Vector{ComplexF32},
    y22::Vector{ComplexF32},
    br::PSY.Transformer3W,
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
    y_shunt = PSY.get_g(br) + im * PSY.get_b(br)
    if primary_available
        primary_ix, star_ix = get_bus_indices(primary_star_arc, num_bus, nr)
        c = 1 / PSY.get_primary_turns_ratio(br)
        Y_t = 1 / (PSY.get_r_primary(br) + PSY.get_x_primary(br) * 1im)
        Y11 = (Y_t * c^2)
        if !isfinite(Y11) || !isfinite(y_shunt)
            error(
                "Data in $(PSY.get_name(br)) is incorrect.
                r_p = $(PSY.get_r_primary(br)), x_p = $(PSY.get_x_primary(br))",
            )
        end

        adj[primary_ix, star_ix] = 1
        adj[star_ix, primary_ix] = -1
        fb[offset_ix + ix + n_entries] = primary_ix
        tb[offset_ix + ix + n_entries] = star_ix

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
        c = 1 / PSY.get_secondary_turns_ratio(br)
        Y_t = 1 / (PSY.get_r_secondary(br) + PSY.get_x_secondary(br) * 1im)
        Y11 = (Y_t * c^2)
        if !isfinite(Y11) || !isfinite(y_shunt)
            error(
                "Data in $(PSY.get_name(br)) is incorrect.
                r_p = $(PSY.get_r_primary(br)), x_p = $(PSY.get_x_primary(br))",
            )
        end

        adj[secondary_ix, star_ix] = 1
        adj[star_ix, secondary_ix] = -1
        fb[offset_ix + ix + n_entries] = secondary_ix
        tb[offset_ix + ix + n_entries] = star_ix

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
        c = 1 / PSY.get_tertiary_turns_ratio(br)
        Y_t = 1 / (PSY.get_r_tertiary(br) + PSY.get_x_tertiary(br) * 1im)
        Y11 = (Y_t * c^2)
        if !isfinite(Y11) || !isfinite(y_shunt)
            error(
                "Data in $(PSY.get_name(br)) is incorrect.
                r_p = $(PSY.get_r_primary(br)), x_p = $(PSY.get_x_primary(br))",
            )
        end

        adj[tertiary_ix, star_ix] = 1
        adj[star_ix, tertiary_ix] = -1
        fb[offset_ix + ix + n_entries] = tertiary_ix
        tb[offset_ix + ix + n_entries] = star_ix

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
    y11::Vector{ComplexF32},
    y12::Vector{ComplexF32},
    y21::Vector{ComplexF32},
    y22::Vector{ComplexF32},
    br::PSY.PhaseShiftingTransformer3W,
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
    y11::Vector{ComplexF32},
    y12::Vector{ComplexF32},
    y21::Vector{ComplexF32},
    y22::Vector{ComplexF32},
    br::PSY.TapTransformer,
    num_bus::Dict{Int, Int},
    branch_ix::Int,
    fb::Vector{Int},
    tb::Vector{Int},
    nr::NetworkReductionData,
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
    y11::Vector{ComplexF32},
    y12::Vector{ComplexF32},
    y21::Vector{ComplexF32},
    y22::Vector{ComplexF32},
    br::PSY.PhaseShiftingTransformer,
    num_bus::Dict{Int, Int},
    branch_ix::Int,
    fb::Vector{Int},
    tb::Vector{Int},
    nr::NetworkReductionData,
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
    Y = PSY.get_impedance_active_power(fa) + im * PSY.get_impedance_reactive_power(fa)
    if !isfinite(Y)
        error(
            "Data in $(PSY.get_name(fa)) is incorrect. Y = $(Y)",
        )
    end
    sb[fa_ix] = bus_no
    ysh[fa_ix] = Y
    return
end

function _ybus!(
    y11::Vector{ComplexF32},
    y12::Vector{ComplexF32},
    y21::Vector{ComplexF32},
    y22::Vector{ComplexF32},
    br::PSY.DiscreteControlledACBranch,
    num_bus::Dict{Int, Int},
    branch_ix::Int,
    fb::Vector{Int},
    tb::Vector{Int},
    nr::NetworkReductionData,
)
    arc = PSY.get_arc(br)
    bus_from_no, bus_to_no = get_bus_indices(arc, num_bus, nr)
    Y_l = (1 / (PSY.get_r(br) + PSY.get_x(br) * 1im))

    if !isfinite(Y_l)
        error(
            "Data in $(PSY.get_name(br)) is incorrect. r = $(PSY.get_r(br)), x = $(PSY.get_x(br))",
        )
    end

    adj[bus_from_no, bus_to_no] = 1
    adj[bus_to_no, bus_from_no] = -1
    fb[branch_ix] = bus_from_no
    tb[branch_ix] = bus_to_no

    y11[branch_ix] = Y_l
    y12[branch_ix] = -Y_l
    y21[branch_ix] = -Y_l
    y22[branch_ix] = Y_l
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
Builds a Ybus from the system. The return is a Ybus Array indexed with the bus numbers and the branch names.

# Arguments
- `check_connectivity::Bool`: Checks connectivity of the network
"""
function Ybus(
    sys::PSY.System;
    check_connectivity::Bool = true,
    make_branch_admittance_matrices::Bool = false,
    network_reductions::Vector{NetworkReduction} = NetworkReduction[],
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

    if check_connectivity && length(bus_lookup) > 1
        subnetworks = assign_reference_buses!(
            find_subnetworks(ybus, axes[1]),
            ref_bus_numbers,
        )
        if length(subnetworks) > 1
            @warn "More than one island found; Network is not connected"
        end
    else
        subnetworks = Dict{Int, Set{Int}}()
    end
    ybus = Ybus(
        ybus,
        adj,
        ref_bus_numbers,
        axes,
        lookup,
        subnetworks,
        nr,
        yft,
        ytf,
        fb,
        tb,
    )

    for nr in network_reductions
        ybus = build_reduced_ybus(ybus, sys, nr)
    end
    return ybus
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
    validate_reduction_type(
        get_reductions(nr_new)[1],
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
    _add_series_branches_to_ybus!(ybus, nr_new.series_branch_map)

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
    elseif !isempty(nr_new.series_branch_map) && !!isempty(nr.series_branch_map)
        error(
            "Cannot compose series branch maps; should not apply multiple reductions that generate series branch maps",
        )
    end
    #TODO - loop through added branches and admittances and modify the Ybus for Ward reduction
    push!(nr.reductions, nr_new.reductions[1])
    bus_ax = setdiff(ybus.axes[1], bus_numbers_to_remove)
    bus_lookup = make_ax_ref(bus_ax)
    bus_ix = [ybus.lookup[1][x] for x in bus_ax]
    adjacency_data = adjacency_data[bus_ix, bus_ix]
    data = data[bus_ix, bus_ix]

    subnetworks = deepcopy(ybus.subnetworks)
    for (k, values) in subnetworks
        if k in bus_numbers_to_remove
            pop!(subnetworks, k)
        else
            for x in values
                if x in bus_numbers_to_remove
                    pop!(subnetworks[k], x)
                end
            end
        end
    end
    #TODO - doesn't yet account for reduction in yft, ytf, fb, tb
    return Ybus(
        data,
        adjacency_data,
        setdiff(ybus.ref_bus_numbers, bus_numbers_to_remove),
        (bus_ax, bus_ax),
        (bus_lookup, bus_lookup),
        subnetworks,
        nr,
        ybus.yft,
        ybus.ytf,
        ybus.fb,
        ybus.tb,
    )
end

function _add_series_branches_to_ybus!(
    ybus::Ybus,
    series_branch_map::Dict{Tuple{Int, Int}, Vector{Any}},
)
    bus_lookup = ybus.lookup[1]
    data = ybus.data
    nrd = get_network_reduction_data(ybus)
    for (equivalent_arc, series_map_entry) in series_branch_map
        series_impedance_from_to = 0.0
        series_impedance_to_from = 0.0
        shunt_admittance = 0.0
        for series_segment in series_map_entry
            segment_series_impedance_from_to,
            segment_series_impedance_to_from,
            segment_shunt_admittance =
                _get_equivalent_line_parameters(ybus, equivalent_arc, series_segment, nrd)
            series_impedance_from_to += segment_series_impedance_from_to
            series_impedance_to_from += segment_series_impedance_to_from
            shunt_admittance += segment_shunt_admittance
        end
        data[bus_lookup[equivalent_arc[1]], bus_lookup[equivalent_arc[2]]] +=
            1 / series_impedance_from_to
        data[bus_lookup[equivalent_arc[2]], bus_lookup[equivalent_arc[1]]] +=
            1 / series_impedance_to_from
        data[bus_lookup[equivalent_arc[1]], bus_lookup[equivalent_arc[1]]] +=
            shunt_admittance / 2
        data[bus_lookup[equivalent_arc[2]], bus_lookup[equivalent_arc[2]]] +=
            shunt_admittance / 2
    end
end

function _get_equivalent_line_parameters(
    ybus::Ybus,
    equivalent_arc::Tuple{Int, Int},
    branch::PSY.Branch,
    network_reduction_data::NetworkReductionData,
)
    fr_bus_no, to_bus_no = get_arc_tuple(branch, network_reduction_data)
    data = ybus.data
    bus_lookup = ybus.lookup[1]
    series_impedance_from_to = 1 / data[bus_lookup[fr_bus_no], bus_lookup[to_bus_no]]
    series_impedance_to_from = 1 / data[bus_lookup[to_bus_no], bus_lookup[fr_bus_no]]
    shunt_admittance = 0.0
    if fr_bus_no ∉ equivalent_arc
        shunt_admittance += data[bus_lookup[fr_bus_no], bus_lookup[fr_bus_no]]
    end
    if to_bus_no ∉ equivalent_arc
        shunt_admittance += data[bus_lookup[to_bus_no], bus_lookup[to_bus_no]]
    end
    return series_impedance_from_to, series_impedance_to_from, shunt_admittance
end

function _get_equivalent_line_parameters(
    ybus::Ybus,
    equivalent_arc::Tuple{Int, Int},
    branches::Set{PSY.Branch},
    network_reduction_data::NetworkReductionData,
)
    error("Implement getting equivalent line parameters for set of parallel branches")
end

function _get_equivalent_line_parameters(
    ybus::Ybus,
    equivalent_arc::Tuple{Int, Int},
    transformer_tuple::Tuple{PSY.ThreeWindingTransformer, Int},
    network_reduction_data::NetworkReductionData,
)
    fr_bus_no, to_bus_no = get_arc_tuple(transformer_tuple, network_reduction_data)
    data = ybus.data
    bus_lookup = ybus.lookup[1]
    series_impedance_from_to = 1 / data[bus_lookup[fr_bus_no], bus_lookup[to_bus_no]]
    series_impedance_to_from = 1 / data[bus_lookup[to_bus_no], bus_lookup[fr_bus_no]]
    shunt_admittance = 0.0
    if fr_bus_no ∉ equivalent_arc
        shunt_admittance += data[bus_lookup[fr_bus_no], bus_lookup[fr_bus_no]]
    end
    if to_bus_no ∉ equivalent_arc
        shunt_admittance += data[bus_lookup[to_bus_no], bus_lookup[to_bus_no]]
    end
    return series_impedance_from_to, series_impedance_to_from, shunt_admittance
end

function _remake_reverse_direct_branch_map!(nr::NetworkReductionData)
    reverse_direct_branch_map = Dict{PSY.Branch, Tuple{Int, Int}}()
    for (k, v) in nr.direct_branch_map
        reverse_direct_branch_map[v] = k
    end
    nr.reverse_direct_branch_map = reverse_direct_branch_map
end
function _remake_reverse_parallel_branch_map!(nr::NetworkReductionData)
    reverse_parallel_branch_map = Dict{PSY.Branch, Tuple{Int, Int}}()
    for (k, v) in nr.parallel_branch_map
        for x in v
            reverse_parallel_branch_map[x] = k
        end
    end
    nr.reverse_parallel_branch_map = reverse_parallel_branch_map
end
function _remake_reverse_series_branch_map!(nr::NetworkReductionData)
    reverse_series_branch_map = Dict{PSY.Branch, Tuple{Int, Int}}()
    for (k, v) in nr.series_branch_map
        for x in v
            reverse_series_branch_map[x] = k
        end
    end
    nr.reverse_series_branch_map = reverse_series_branch_map
end
function _remake_reverse_transformer3W_map!(nr::NetworkReductionData)
    reverse_transformer3W_map =
        Dict{Tuple{PSY.ThreeWindingTransformer, Int}, Tuple{Int, Int}}()
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
