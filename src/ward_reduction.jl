"""
Builds a NetworkReduction corresponding to Ward reduction

# Arguments
- `sys::System`
- `study_buses::Vector{Int}`: Bus numbers corresponding to the area of study (the retained area)
"""
function get_ward_reduction(
    sys::PSY.System,
    study_buses::Vector{Int};
    prior_reduction::NetworkReduction = NetworkReduction(),
)
    validate_reduction_type(NetworkReductionTypes.WARD, get_reduction_type(prior_reduction))
    _validate_study_buses(sys, study_buses, prior_reduction)
    y_bus = Ybus(sys; network_reduction = prior_reduction, check_connectivity = false)
    Z_full = KLU.solve!(klu(y_bus.data), Matrix(one(y_bus.data)))       #TODO: change implementation for large systems (row by row)
    boundary_buses = Vector{Int64}()
    for branch in get_ac_branches(sys, prior_reduction.removed_branches)
        (from_bus, to_bus) = get_arc_tuple(branch)
        if (from_bus ∈ study_buses) && (to_bus ∉ study_buses)
            push!(boundary_buses, from_bus)
        end
        if (to_bus ∈ study_buses) && (from_bus ∉ study_buses)
            push!(boundary_buses, to_bus)
        end
    end
    all_buses =
        [PSY.get_number(b) for b in get_buses(sys, prior_reduction.bus_reduction_map)]
    external_buses = setdiff(all_buses, study_buses)
    boundary_buses = unique(boundary_buses)
    n_external = length(external_buses)
    n_boundary = length(boundary_buses)

    retained_branches = Set{String}()
    removed_branches = Set{String}()
    for branch in get_ac_branches(sys, prior_reduction.removed_branches)
        (from_bus, to_bus) = get_arc_tuple(branch)
        if (from_bus ∈ external_buses) && (to_bus ∈ external_buses) ||
           (from_bus ∈ boundary_buses) && (to_bus ∈ external_buses) ||
           (from_bus ∈ external_buses) && (to_bus ∈ boundary_buses)
            push!(removed_branches, PSY.get_name(branch))
        else
            push!(retained_branches, PSY.get_name(branch))
        end
    end

    bus_reduction_map_index = Dict{Int, Set{Int}}(k => Set{Int}() for k in study_buses)
    bus_lookup = y_bus.lookup[1]    #y_bus and Z have same lookup
    for b in external_buses
        boundary_bus_indices = [bus_lookup[x] for x in boundary_buses]
        boundary_bus_numbers = [x for x in boundary_buses]
        row_index = bus_lookup[b]
        Z_row_boundary = abs.(Z_full[row_index, boundary_bus_indices])
        closest_boundary_bus = boundary_bus_numbers[argmin(Z_row_boundary)]
        push!(bus_reduction_map_index[closest_boundary_bus], b)
    end
    reverse_bus_search_map =
        _make_reverse_bus_search_map(bus_reduction_map_index, length(all_buses))

    #Populate matrices for computing external equivalent
    y_ee = SparseArrays.spzeros(ComplexF64, n_external, n_external)
    for (ix, i) in enumerate(external_buses)
        for (jx, j) in enumerate(external_buses)
            y_ee[ix, jx] = y_bus[i, j]
        end
    end
    y_be = SparseArrays.spzeros(ComplexF64, n_boundary, n_external)
    for (ix, i) in enumerate(boundary_buses)
        for (jx, j) in enumerate(external_buses)
            y_be[ix, jx] = y_bus[i, j]
        end
    end
    y_eb = SparseArrays.spzeros(ComplexF64, n_external, n_boundary)
    for (ix, i) in enumerate(external_buses)
        for (jx, j) in enumerate(boundary_buses)
            y_eb[ix, jx] = y_bus[i, j]
        end
    end

    # Eq. (2.16) from  https://core.ac.uk/download/pdf/79564835.pdf
    y_eq = y_be * KLU.solve!(klu(y_ee), Matrix(y_eb))

    added_branches = Vector{PSY.ACTransmission}()
    added_admittances = Vector{PSY.FixedAdmittance}()
    virtual_admittance_name_index = 1
    virtual_branch_name_index = 1
    #Loop upper diagonal of Yeq
    for ix in 1:length(boundary_buses)
        for jx in ix:length(boundary_buses)
            if y_eq[ix, jx] != 0.0
                if ix == jx
                    bus = collect(
                        PSY.get_components(
                            x -> PSY.get_number(x) == boundary_buses[ix],
                            PSY.ACBus,
                            sys,
                        ),
                    )[1]
                    virtual_admittance = PSY.FixedAdmittance(;
                        name = "virtual_admittance_$(virtual_admittance_name_index)",
                        available = true,
                        bus = bus,
                        Y = y_eq[ix, jx],
                    )
                    push!(added_admittances, virtual_admittance)
                    virtual_admittance_name_index += 1
                else
                    to_bus = collect(
                        PSY.get_components(
                            x -> PSY.get_number(x) == boundary_buses[ix],
                            PSY.ACBus,
                            sys,
                        ),
                    )[1]
                    from_bus = collect(
                        PSY.get_components(
                            x -> PSY.get_number(x) == boundary_buses[jx],
                            PSY.ACBus,
                            sys,
                        ),
                    )[1]
                    virtual_branch = PSY.Line(;
                        name = "virtual_branch_$(virtual_branch_name_index)",
                        available = true,
                        active_power_flow = 0.0,
                        reactive_power_flow = 0.0,
                        arc = PSY.Arc(; from = from_bus, to = to_bus),
                        r = -1 * real(y_eq[ix, jx]),
                        x = -1 * imag(y_eq[ix, jx]),
                        b = (from = 0.0, to = 0.0),
                        rating = 100.0,
                        angle_limits = (min = (-pi / 3), max = (pi / 3)),
                        g = (from = 0.0, to = 0.0),
                    )
                    push!(added_branches, virtual_branch)
                    virtual_branch_name_index += 1
                end
            end
        end
    end
    new_reduction = NetworkReduction(;
        bus_reduction_map = bus_reduction_map_index,
        reverse_bus_search_map = reverse_bus_search_map,
        removed_branches = removed_branches,
        retained_branches = retained_branches,
        added_branches = added_branches,
        added_admittances = added_admittances,
        reduction_type = [NetworkReductionTypes.WARD],
    )
    if isempty(prior_reduction)
        return new_reduction
    else
        return compose_reductions(prior_reduction, new_reduction, length(all_buses))
    end
end

function _validate_study_buses(
    sys::PSY.System,
    study_buses::Vector{Int},
    network_reduction::NetworkReduction,
)
    buses = get_buses(sys, network_reduction.bus_reduction_map)
    branches = get_ac_branches(sys, network_reduction.removed_branches)

    bus_numbers = PSY.get_number.(buses)
    for b in study_buses
        b ∉ bus_numbers && throw(IS.DataFormatError("Study bus $b not found in system"))
    end
    M, _ = calculate_adjacency(branches, buses, network_reduction)
    bus_ax = PSY.get_number.(buses)
    sub_nets = find_subnetworks(M, bus_ax)
    if length(sub_nets) > 1
        @warn "System contains multiple islands"
    end

    slack_bus_numbers =
        [PSY.get_number(n) for n in buses if PSY.get_bustype(n) == ACBusTypes.REF]
    for (_, v) in sub_nets
        all_in = all(x -> x in Set(v), study_buses)
        none_in = all(x -> !(x in Set(v)), study_buses)
        if all_in
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
            if sb in v && sb ∉ study_buses && !(none_in)
                throw(
                    IS.DataFormatError(
                        "Slack bus $sb must be included in the study buses for an area that is partially reduced",
                    ),
                )
            end
        end
    end
end
