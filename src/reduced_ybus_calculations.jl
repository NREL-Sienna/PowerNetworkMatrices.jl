struct ReducedYbus{Ax, L <: NTuple{2, Dict}} <: PowerNetworkMatrix{ComplexF64}
    data::SparseArrays.SparseMatrixCSC{ComplexF64, Int}
    axes::Ax
    lookup::L
end

function WardReduction(y_bus::Ybus, sys::PSY.System, study_buses::Vector{Int64})
    reduced_axes = (study_buses, study_buses)
    reduced_bus_lookup = make_ax_ref(study_buses)
    reduced_lookup = (reduced_bus_lookup, reduced_bus_lookup)
    reduced_bus_count = length(study_buses)
    reduced_data = SparseArrays.spzeros(ComplexF64, reduced_bus_count, reduced_bus_count)
    Ybus_reduced = ReducedYbus(reduced_data, reduced_axes, reduced_lookup)

    #Match original entries for internal and boundary buses.
    for i in study_buses
        for j in study_buses
            Ybus_reduced[i, j] = y_bus[i, j]
        end
    end

    #Find boundary buses (subset of study buses)
    boundary_buses = Int64[]
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
    n_boundary = length(boundary_buses)
    n_external = length(external_buses)

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

    # Modify boundary entries with equivalent 
    for (ix, i) in enumerate(boundary_buses)
        for (jx, j) in enumerate(boundary_buses)
            Ybus_reduced[i, j] -= y_eq[ix, jx]
        end
    end

    return Ybus_reduced
end
