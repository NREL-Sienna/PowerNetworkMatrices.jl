@testset "Test 3W Integration to Ybus Matrix" begin
    sys = System("test_data/data/6bus_3w_system_updated.raw")
    Ybus_pnm = Ybus(sys)

    data = readlines("test_data/data/6bus_3w_system_updated_Ymx.txt")
    from_buses, to_buses, y_values = Int[], Int[], ComplexF64[]
    for line in data
        vals = split(strip(line), r",\s*")
        from_bus = parse(Int, vals[1])
        to_bus = parse(Int, vals[2])
        y_real = parse(Float64, vals[3])
        y_imag = parse(Float64, vals[4])
        push!(from_buses, from_bus)
        push!(to_buses, to_bus)
        push!(y_values, y_real + im * y_imag)
    end
    unique_buses = sort(unique(vcat(from_buses, to_buses)))
    bus_index = Dict(bus => i for (i, bus) in enumerate(unique_buses))
    row_indices = [bus_index[f] for f in from_buses]
    col_indices = [bus_index[t] for t in to_buses]
    println("PSSE Complex Y-bus matrix:")
    Ybus_psse = sparse(row_indices, col_indices, y_values)
    # Insert missing PSSE Ymatrix elements
    Ybus_psse[1,1] = 8.43881+0.0535936im
    Ybus_psse[1,2] = -8.43881+0.00712136im

    nrows = size(Ybus_pnm)[1]
    ncols = size(Ybus_pnm)[2]

    for i in range(1, nrows)
        for j in range(1, ncols)
            @test isapprox(Ybus_pnm.data[i, j], Ybus_psse[i, j], atol = 1e-4)
        end
    end
end
