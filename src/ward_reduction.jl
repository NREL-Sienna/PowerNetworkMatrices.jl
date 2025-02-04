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

    virtual_admittances = Vector{Tuple{Int64, Int64, ComplexF64}}()
    for (ix, x) in enumerate(boundary_buses)
        for (iy, y) in enumerate(boundary_buses)
            push!(virtual_admittances, (x, y, y_eq[ix, iy]))
        end
    end

    return NetworkReduction(;
        bus_reduction_map = bus_reduction_map_index,
        reverse_bus_search_map = reverse_bus_search_map,
        removed_branches = removed_branches,
        retained_branches = retained_branches,
        virtual_admittances = virtual_admittances,
    )
end
