"""
Nodal admittance matrix (Ybus) is an N x N matrix describing a power system with N buses. It represents the nodal admittance of the buses in a power system.

The Ybus Struct is indexed using the Bus Numbers, no need for them to be sequential

The fields yft and ytf are the branch admittance matrices for the from-to and to-from branch admittances respectively. The rows correspond to branches and the columns to buses.
The matrix columns are mapped to buses using fb, tb arrays of the matrix columns that correspond to the `from` and `to` buses. 
Using yft, ytf, and the voltage vector, the branch currents and power flows can be calculated.
"""
struct Ybus{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{ComplexF64}
    data::SparseArrays.SparseMatrixCSC{ComplexF64, Int}
    axes::Ax
    lookup::L
    yft::Union{SparseArrays.SparseMatrixCSC{ComplexF64, Int}, Nothing}
    ytf::Union{SparseArrays.SparseMatrixCSC{ComplexF64, Int}, Nothing}
    fb::Union{Vector{Int64}, Nothing}
    tb::Union{Vector{Int64}, Nothing}
end

function _ybus!(
    y11::Vector{ComplexF64},
    y12::Vector{ComplexF64},
    y21::Vector{ComplexF64},
    y22::Vector{ComplexF64},
    br::PSY.ACBranch,
    num_bus::Dict{Int, Int},
    branch_ix::Int64,
    fb::Vector{Int64},
    tb::Vector{Int64},
)
    arc = PSY.get_arc(br)
    bus_from_no = num_bus[arc.from.number]
    bus_to_no = num_bus[arc.to.number]
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
    br::PSY.DynamicBranch,
    num_bus::Dict{Int, Int},
    branch_ix::Int64,
    fb::Vector{Int64},
    tb::Vector{Int64},
)
    _ybus!(y11, y12, y21, y22, br.branch, num_bus, branch_ix, fb, tb)
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
)
    arc = PSY.get_arc(br)
    bus_from_no = num_bus[arc.from.number]
    bus_to_no = num_bus[arc.to.number]
    fb[branch_ix] = bus_from_no
    tb[branch_ix] = bus_to_no
    Y_t = 1 / (PSY.get_r(br) + PSY.get_x(br) * 1im)
    Y11 = Y_t
    b = PSY.get_primary_shunt(br)
    if !isfinite(Y11) || !isfinite(Y_t) || !isfinite(b)
        error(
            "Data in $(PSY.get_name(br)) is incorrect. r = $(PSY.get_r(br)), x = $(PSY.get_x(br))",
        )
    end

    y11[branch_ix] = Y11 - (1im * b)
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
    branch_ix::Int64,
    fb::Vector{Int64},
    tb::Vector{Int64},
    stb::Int64,
)
    primary_secondary_arc = PSY.get_primary_secondary_arc(br)
    secondary_tertiary_arc = PSY.get_secondary_tertiary_arc(br)
    primary_tertiary_arc = PSY.get_primary_tertiary_arc(br)
    start_bus = PSY.get_star_bus(br)

    bus_ps_from_no = num_bus[primary_secondary_arc.from.number]
    bus_st_from_no = num_bus[secondary_tertiary_arc.from.number]
    bus_pt_to_no = num_bus[primary_tertiary_arc.to.number]

    fb[branch_ix + stb] = bus_ps_from_no
    fb[branch_ix + stb + 1] = bus_st_from_no
    fb[branch_ix + stb + 2] = bus_pt_to_no

    tb[branch_ix + stb] = num_bus[start_bus.number]
    tb[branch_ix + stb + 1] = num_bus[start_bus.number]
    tb[branch_ix + stb + 2] = num_bus[start_bus.number]

    Y_t1 = 1 / (PSY.get_r_primary(br) + PSY.get_x_primary(br) * 1im)
    Y11 = Y_t1
    Y_t2 = 1 / (PSY.get_r_secondary(br) + PSY.get_x_secondary(br) * 1im)
    Y22 = Y_t2
    Y_t3 = 1 / (PSY.get_r_tertiary(br) + PSY.get_x_tertiary(br) * 1im)
    Y33 = Y_t3

    b = PSY.get_b(br) # shunt susceptance (star bus to ground)

    if !isfinite(Y11) || !isfinite(Y22) || !isfinite(Y33) || !isfinite(b)
        error(
            "Data in $(PSY.get_name(br)) is incorrect. 
            r_p = $(PSY.get_r_primary(br)), x_p = $(PSY.get_x_primary(br)), 
            r_s = $(PSY.get_r_secondary(br)), x_s = $(PSY.get_x_secondary(br)), 
            r_t = $(PSY.get_r_tertiary(br)), x_t = $(PSY.get_x_tertiary(br))",
        )
    end

    y11[branch_ix + stb] = Y11 - (1im * b)
    y11[branch_ix + stb + 1] = Y22 - (1im * b)
    y11[branch_ix + stb + 2] = Y33 - (1im * b)

    y12[branch_ix + stb] = -Y_t1
    y12[branch_ix + stb + 1] = -Y_t2
    y12[branch_ix + stb + 2] = -Y_t3

    y21[branch_ix + stb] = -Y_t1
    y21[branch_ix + stb + 1] = -Y_t2
    y21[branch_ix + stb + 2] = -Y_t3

    y22[branch_ix + stb] = Y11
    y22[branch_ix + stb + 1] = Y22
    y22[branch_ix + stb + 2] = Y33

    return
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
)
    arc = PSY.get_arc(br)
    bus_from_no = num_bus[arc.from.number]
    bus_to_no = num_bus[arc.to.number]
    fb[branch_ix] = bus_from_no
    tb[branch_ix] = bus_to_no

    Y_t = 1 / (PSY.get_r(br) + PSY.get_x(br) * 1im)
    c = 1 / PSY.get_tap(br)
    b = PSY.get_primary_shunt(br)

    Y11 = (Y_t * c^2)
    y11[branch_ix] = Y11
    Y12 = (-Y_t * c)
    if !isfinite(Y11) || !isfinite(Y12) || !isfinite(b)
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
)
    arc = PSY.get_arc(br)
    bus_from_no = num_bus[arc.from.number]
    bus_to_no = num_bus[arc.to.number]
    fb[branch_ix] = bus_from_no
    tb[branch_ix] = bus_to_no
    Y_t = 1 / (PSY.get_r(br) + PSY.get_x(br) * 1im)
    tap = (PSY.get_tap(br) * exp(PSY.get_α(br) * 1im))
    c_tap = (PSY.get_tap(br) * exp(-1 * PSY.get_α(br) * 1im))
    b = PSY.get_primary_shunt(br)
    Y11 = (Y_t / abs(tap)^2)
    if !isfinite(Y11) || !isfinite(Y_t) || !isfinite(b * c_tap)
        error(
            "Data in $(PSY.get_name(br)) is incorrect. r = $(PSY.get_r(br)), x = $(PSY.get_x(br))",
        )
    end

    y11[branch_ix] = Y11 - (1im * b)
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
    fa::Union{PSY.FixedAdmittance, PSY.SwitchedAdmittance},
    num_bus::Dict{Int, Int},
    fa_ix::Int64,
    sb::Vector{Int64},
)
    bus = PSY.get_bus(fa)
    bus_no = num_bus[PSY.get_number(bus)]
    sb[fa_ix] = bus_no
    if !isfinite(fa.Y)
        error(
            "Data in $(PSY.get_name(fa)) is incorrect. Y = $(fa.Y)",
        )
    end
    ysh[fa_ix] = fa.Y
    return
end

function _buildybus(
    branches,
    # transformer_3W,
    buses::Vector{PSY.ACBus},
    fixed_admittances::Vector{PSY.FixedAdmittance},
    switched_admittances::Vector{PSY.SwitchedAdmittance},
)
    num_bus = Dict{Int, Int}()

    branchcount = length(branches)
    fa_count = length(fixed_admittances)
    sa_count = length(switched_admittances)
    fb = zeros(Int64, branchcount)
    tb = zeros(Int64, branchcount)
    sb = zeros(Int64, fa_count + sa_count)

    for (ix, b) in enumerate(buses)
        num_bus[PSY.get_number(b)] = ix
    end

    y11 = zeros(ComplexF64, branchcount)
    y12 = zeros(ComplexF64, branchcount)
    y21 = zeros(ComplexF64, branchcount)
    y22 = zeros(ComplexF64, branchcount)
    ysh = zeros(ComplexF64, fa_count + sa_count)

    ix = 0
    for b in branches
        if PSY.get_name(b) == "init"
            throw(DataFormatError("The data in Branch is invalid"))
        end
        ix = _next_branch_number!(b, ix)
        PSY.get_available(b) && _ybus!(y11, y12, y21, y22, b, num_bus, ix, fb, tb)
    end

    # stb = 0
    # for (ix, b) in enumerate(transformer_3W)
    #     if PSY.get_name(b) == "init"
    #         throw(DataFormatError("The data in Transformer3W is invalid"))
    #     end
    #     PSY.get_available(b) &&
    #         _ybus!(y11, y12, y21, y22, b, num_bus, ix + branchcount_no_3w, fb, tb, stb)

    #     stb = stb + 2
    # end

    for (ix, fa) in enumerate([fixed_admittances; switched_admittances])
        PSY.get_available(fa) && _ybus!(ysh, fa, num_bus, ix, sb)
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
Builds a Ybus from a collection of buses and branches. The return is a Ybus Array indexed with the bus numbers and the branch names.

# Arguments
- `check_connectivity::Bool`: Checks connectivity of the network using Depth First Search (DFS)
"""
function Ybus(
    branches::Vector,
    buses::Vector{PSY.ACBus},
    # transformer_3W::Vector{PSY.Transformer3W} = Vector{PSY.Transformer3W}(),
    fixed_admittances::Vector{PSY.FixedAdmittance} = Vector{PSY.FixedAdmittance}(),
    switched_admittances::Vector{PSY.SwitchedAdmittance} = Vector{PSY.SwitchedAdmittance}();
    check_connectivity::Bool = true,
    make_branch_admittance_matrices::Bool = false,
)
    bus_ax = PSY.get_number.(buses)
    axes = (bus_ax, bus_ax)
    bus_lookup = make_ax_ref(bus_ax)
    busnumber = length(buses)
    look_up = (bus_lookup, bus_lookup)
    y11, y12, y21, y22, ysh, fb, tb, sb =
    # _buildybus(branches, transformer_3W, buses, fixed_admittances, switched_admittances)
        _buildybus(branches, buses, fixed_admittances, switched_admittances)
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
    return Ybus(ybus, axes, look_up, yft, ytf, fb, tb)
end

"""
Builds a Ybus from the system. The return is a Ybus Array indexed with the bus numbers and the branch names.

# Arguments
- `check_connectivity::Bool`: Checks connectivity of the network
"""
function Ybus(sys::PSY.System; kwargs...)
    branches = get_ac_branches(sys)
    # transformer_3W = get_transformers_3w(sys)
    buses = get_buses(sys)
    fixed_admittances = collect(PSY.get_components(PSY.FixedAdmittance, sys))
    switched_admittances = collect(PSY.get_components(PSY.SwitchedAdmittance, sys))
    return Ybus(
        branches,
        buses,
        # transformer_3W,
        fixed_admittances,
        switched_admittances;
        kwargs...,
    )
end

"""
Nodal incidence matrix (Adjacency) is an N x N matrix describing a power system with N buses. It represents the directed connectivity of the buses in a power system.

The Adjacency Struct is indexed using the Bus Numbers, no need for them to be sequential
"""
struct Adjacency{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{Int}
    data::SparseArrays.SparseMatrixCSC{Int, Int}
    axes::Ax
    lookup::L
end

"""
Builds a Adjacency from a collection of buses and branches. The return is an N x N Adjacency Array indexed with the bus numbers.

# Keyword arguments
- `check_connectivity::Bool`: Checks connectivity of the network using Goderya's algorithm
"""
function Adjacency(branches, nodes; check_connectivity::Bool = true, kwargs...)
    buscount = length(nodes)
    bus_ax = PSY.get_number.(nodes)
    axes = (bus_ax, bus_ax)
    bus_lookup = make_ax_ref(bus_ax)
    look_up = (bus_lookup, bus_lookup)

    a = SparseArrays.spzeros(Int, buscount, buscount)

    for b in branches
        (fr_b, to_b) = get_bus_indices(b, bus_lookup)
        a[fr_b, to_b] = 1
        a[to_b, fr_b] = -1
        a[fr_b, fr_b] = 1
        a[to_b, to_b] = 1
    end

    if check_connectivity
        connected = validate_connectivity(a, nodes, bus_lookup; kwargs...)
        !connected && throw(DataFormatError("Network not connected"))
    end

    return Adjacency(a, axes, look_up)
end

"""
Builds a Adjacency from the system. The return is an N x N Adjacency Array indexed with the bus numbers.

# Keyword arguments
- `check_connectivity::Bool`: Checks connectivity of the network using Goderya's algorithm
- `connectivity_method::Function = goderya_connectivity`: method (`goderya_connectivity` or `dfs_connectivity`) for connectivity validation
"""
function Adjacency(sys::PSY.System; check_connectivity::Bool = true, kwargs...)
    nodes = sort!(
        collect(
            PSY.get_components(
                x -> PSY.get_bustype(x) != ACBusTypes.ISOLATED,
                PSY.ACBus,
                sys,
            ),
        );
        by = x -> PSY.get_number(x),
    )
    branches = PSY.get_components(PSY.get_available, PSY.Branch, sys)
    return Adjacency(branches, nodes; check_connectivity = check_connectivity, kwargs...)
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
Finds the set of bus numbers that belong to each connnected component in the System
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
