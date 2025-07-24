struct WardReduction <: NetworkReduction
    study_buses::Vector{Int}
end
get_study_buses(nr::WardReduction) = nr.study_buses

"""
Builds a NetworkReduction corresponding to Ward reduction

# Arguments
- `sys::System`
- `study_buses::Vector{Int}`: Bus numbers corresponding to the area of study (the retained area)
"""
function get_ward_reduction(
    y_bus::Ybus,
    study_buses::Vector{Int},
    reduction::WardReduction,
)
    Z_full = KLU.solve!(klu(y_bus.data), Matrix{ComplexF64}(one(y_bus.data)))       #TODO: change implementation for large systems (row by row)
    A = IncidenceMatrix(y_bus)
    boundary_buses = Set{Int}()
    removed_arcs = Set{Tuple{Int, Int}}()
    for arc in A.axes[1]
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
    all_buses = y_bus.axes[1]
    external_buses = setdiff(all_buses, study_buses)
    boundary_buses = unique(boundary_buses)
    n_external = length(external_buses)
    n_boundary = length(boundary_buses)

    bus_reduction_map_index = Dict{Int, Set{Int}}(k => Set{Int}() for k in study_buses)
    bus_lookup = y_bus.lookup[1]    #y_bus and Z have same lookup

    added_branch_map = Dict{Tuple{Int, Int}, PSY.Line}()
    added_admittance_map = Dict{Int, PSY.FixedAdmittance}()
    if isempty(boundary_buses)
        first_ref_study_bus = findfirst(x -> x ∈ y_bus.ref_bus_numbers, study_buses)
        @error "no boundary buses found; cannot make bus_reduction_map based on impedance based criteria. mapping all external buses to the first reference bus ($first_ref_study_bus)"
        bus_reduction_map_index[first_ref_study_bus] = Set(external_buses)
    else
        for b in external_buses
            boundary_bus_indices = [bus_lookup[x] for x in boundary_buses]
            boundary_bus_numbers = [x for x in boundary_buses]
            row_index = bus_lookup[b]
            Z_row_boundary = abs.(Z_full[row_index, boundary_bus_indices])
            closest_boundary_bus = boundary_bus_numbers[argmin(Z_row_boundary)]
            push!(bus_reduction_map_index[closest_boundary_bus], b)
        end
    end
    reverse_bus_search_map =
        _make_reverse_bus_search_map(bus_reduction_map_index, length(all_buses))

    #Populate matrices for computing external equivalent
    y_ee = SparseArrays.spzeros(ComplexF32, n_external, n_external)
    for (ix, i) in enumerate(external_buses)
        for (jx, j) in enumerate(external_buses)
            y_ee[ix, jx] = y_bus[i, j]
        end
    end
    y_be = SparseArrays.spzeros(ComplexF32, n_boundary, n_external)
    for (ix, i) in enumerate(boundary_buses)
        for (jx, j) in enumerate(external_buses)
            y_be[ix, jx] = y_bus[i, j]
        end
    end
    y_eb = SparseArrays.spzeros(ComplexF32, n_external, n_boundary)
    for (ix, i) in enumerate(external_buses)
        for (jx, j) in enumerate(boundary_buses)
            y_eb[ix, jx] = y_bus[i, j]
        end
    end

    # Eq. (2.16) from  https://core.ac.uk/download/pdf/79564835.pdf
    y_eq = y_be * KLU.solve!(klu(y_ee), Matrix{Complex{Float64}}(y_eb))

    virtual_admittance_name_index = 1
    virtual_branch_name_index = 1
    #Loop upper diagonal of Yeq
    for ix in 1:length(boundary_buses)
        for jx in ix:length(boundary_buses)
            bus_ix = boundary_buses[ix]
            bus_jx = boundary_buses[jx]
            if y_eq[ix, jx] != 0.0
                if ix == jx
                    virtual_admittance = PSY.FixedAdmittance(;
                        name = "virtual_admittance_$(virtual_admittance_name_index)",
                        available = true,
                        bus = PSY.ACBus(nothing),
                        Y = y_eq[ix, jx],
                    )
                    added_admittance_map[ix] = virtual_admittance
                    virtual_admittance_name_index += 1
                else
                    #check if the arc of virtual line is already existing so we don't add an additional arc
                    if (bus_ix, bus_jx) ∈ A.axes[1]
                        arc_key = (bus_ix, bus_jx)
                    else
                        arc_key = (bus_jx, bus_ix)
                    end
                    virtual_branch = PSY.Line(;
                        name = "virtual_branch_$(virtual_branch_name_index)",
                        available = true,
                        active_power_flow = 0.0,
                        reactive_power_flow = 0.0,
                        arc = PSY.Arc(; from = PSY.ACBus(nothing), to = PSY.ACBus(nothing)),
                        r = -1 * real(y_eq[ix, jx]),
                        x = -1 * imag(y_eq[ix, jx]),
                        b = (from = 0.0, to = 0.0),
                        rating = 100.0,
                        angle_limits = (min = (-pi / 3), max = (pi / 3)),
                        g = (from = 0.0, to = 0.0),
                    )
                    added_branch_map[arc_key] = virtual_branch
                    virtual_branch_name_index += 1
                end
            end
        end
    end
    return NetworkReductionData(;
        bus_reduction_map = bus_reduction_map_index,
        reverse_bus_search_map = reverse_bus_search_map,
        removed_arcs = removed_arcs,
        added_branch_map = added_branch_map,
        added_admittance_map = added_admittance_map,
        reductions = NetworkReduction[reduction],
    )
end
