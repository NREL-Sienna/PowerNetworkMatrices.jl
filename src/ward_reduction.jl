"""
Builds a NetworkReduction corresponding to Ward reduction

# Arguments
- `sys::System`
- `study_buses::Vector{Int}`: Bus numbers corresponding to the area of study (the retained area)
"""
function get_ward_reduction(sys::PSY.System, study_buses::Vector{Int})
    boundary_buses = Vector{Int64}()
    for branch in PSY.get_components(PSY.ACBranch, sys)
        from_bus = PSY.get_number(PSY.get_from(PSY.get_arc(branch)))
        to_bus = PSY.get_number(PSY.get_to(PSY.get_arc(branch)))
        if (from_bus ∈ study_buses) && (to_bus ∉ study_buses)
            push!(boundary_buses, from_bus)
        end
        if (to_bus ∈ study_buses) && (from_bus ∉ study_buses)
            push!(boundary_buses, to_bus)
        end
    end
    all_buses = [PSY.get_number(b) for b in PSY.get_components(PSY.ACBus, sys)]
    external_buses = setdiff(all_buses, study_buses)
    boundary_buses = unique(boundary_buses)
    n_external = length(external_buses)
    n_boundary = length(boundary_buses)

    retained_branches = Set{String}()
    removed_branches = Set{String}()
    for branch in PSY.get_components(PSY.Branch, sys)
        arc = PSY.get_arc(branch)
        from_bus = PSY.get_number(PSY.get_from(arc))
        to_bus = PSY.get_number(PSY.get_to(arc))
        if (from_bus ∈ external_buses) && (to_bus ∈ external_buses) ||
           (from_bus ∈ boundary_buses) && (to_bus ∈ external_buses) ||
           (from_bus ∈ external_buses) && (to_bus ∈ boundary_buses)
            push!(removed_branches, PSY.get_name(branch))
        else
            push!(retained_branches, PSY.get_name(branch))
        end
    end
    y_bus = Ybus(sys)
    Z_full = KLU.solve!(klu(y_bus.data), Matrix(one(y_bus.data))) #TODO: change implementation for large systems (row by row)
    @assert length(all_buses) < 1000

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

    added_branches = Vector{PSY.Branch}()
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

    return NetworkReduction(;
        bus_reduction_map = bus_reduction_map_index,
        reverse_bus_search_map = reverse_bus_search_map,
        removed_branches = removed_branches,
        retained_branches = retained_branches,
        added_branches = added_branches,
        added_admittances = added_admittances,
    )
end
