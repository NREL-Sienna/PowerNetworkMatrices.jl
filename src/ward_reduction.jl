struct WardReduction <: NetworkReduction
    study_buses::Vector{Int}
end
get_study_buses(nr::WardReduction) = nr.study_buses

function get_reduction(
    ybus::Ybus,
    ::PSY.System,
    reduction::WardReduction,
)
    study_buses = get_study_buses(reduction)
    return get_ward_reduction(ybus, study_buses, reduction)
end

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
    _validate_study_buses(y_bus, study_buses)
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

    added_branch_map = Dict{Tuple{Int, Int}, PSY.Line}()
    added_admittance_map = Dict{Int, PSY.FixedAdmittance}()
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
    slack_bus_numbers = ybus.ref_bus_numbers
    if isempty(ybus.subnetworks)
        @warn "Skipping additional data checks because subnetworks are not computed."
    else
        sub_networks = ybus.subnetworks
        for (_, v) in sub_networks
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
end
