using SparseArrays

function parse_psse_ybus(path)
    data = readlines(path)
    row_buses, col_buses, y_values = Int[], Int[], ComplexF64[]
    reduced_bus_pairs = Vector{Vector{Int}}()
    values_done = false
    for line in data
        if line == " ZERO IMPEDANCE LINE CONNECTED BUSES (BUS IN MATRIX LISTED FIRST):"
            values_done = true
        end
        if !isempty(line) &&
           line != " ZERO IMPEDANCE LINE CONNECTED BUSES (BUS IN MATRIX LISTED FIRST):"
            if values_done == true
                vals = split(strip(line), r",\s*")
                push!(reduced_bus_pairs, [parse(Int, x) for x in vals])
            else
                vals = split(strip(line), r",\s*")
                row_bus = parse(Int, vals[1])
                col_bus = parse(Int, vals[2])
                y_real = parse(Float64, vals[3])
                y_imag = parse(Float64, vals[4])
                push!(row_buses, row_bus)
                push!(col_buses, col_bus)
                push!(y_values, y_real + im * y_imag)
            end
        end
    end
    unique_buses = sort(unique(vcat(row_buses, col_buses)))
    bus_to_ix_map = Dict(bus => i for (i, bus) in enumerate(unique_buses))
    row_indices = [bus_to_ix_map[f] for f in row_buses]
    col_indices = [bus_to_ix_map[t] for t in col_buses]
    return sparse(row_indices, col_indices, y_values),
    bus_to_ix_map,
    row_buses,
    col_buses,
    y_values,
    reduced_bus_pairs
end

function _test_psse_reduction_row(psse_row, reverse_bus_search_map)
    reduced = zeros(length(psse_row))
    reduced_to = zeros(length(psse_row))
    for (ix, x) in enumerate(psse_row)
        if haskey(reverse_bus_search_map, x)
            reduced[ix] = 1
            reduced_to[ix] = reverse_bus_search_map[x]
        end
    end
    @test sum(reduced) == length(psse_row) - 1   #all but one entry should be reduced 
    ix = findfirst(isequal(0), reduced)
    splice!(reduced_to, ix)
    @test all(reduced_to .== psse_row[ix])       #all entries should be reduced to the one entry that was not reduced
end

@testset "14 bus system" begin
    sys = System(joinpath(pwd(), "test", "test_data", "modified_14bus_system.raw"))
    Ybus_pnm = Ybus(sys; check_connectivity = false)
    ref_bus_numbers = [
        get_number(x) for
        x in get_components(x -> get_bustype(x) == PSY.ACBusTypes.REF, ACBus, sys)
    ]
    skip_indices = indexin(ref_bus_numbers, Ybus_pnm.axes[1])
    n_ref_bus_elements = nnz(Ybus_pnm.data[skip_indices, :])
    nr = Ybus_pnm.network_reduction
    # Test the NetworkReduction:
    #Test breaker/switch reduction: 
    @test length(nr.bus_reduction_map) == 14
    @test length(nr.reverse_bus_search_map) == 2
    #Test direct branch maps (1:1):
    @test length(keys(nr.direct_branch_map)) ==
          length(keys(nr.reverse_direct_branch_map)) == 13
    #Test double circuit maps:
    @test length(nr.parallel_branch_map) == 1
    @test length(nr.parallel_branch_map[(106, 111)]) == 2
    @test length(nr.reverse_parallel_branch_map) == 2
    #Test 3WT maps (1:1):
    @test length(nr.transformer3W_map) == length(nr.reverse_transformer3W_map) == 6
    #Test adjacency and Ybus have same non-zero elements
    findnz(Ybus_pnm.data)[1] == findnz(Ybus_pnm.adjacency_data)[1]
    findnz(Ybus_pnm.data)[2] == findnz(Ybus_pnm.adjacency_data)[2]

    # Compare with PSSE
    Ybus_psse, b_ix_psse, row_buses, col_buses, y_values, reduced_bus_pairs_psse =
        parse_psse_ybus("test/test_data/modified_14bus_system_Ymx_SWTBRK.txt")  # NOTE: manually changed the ybus numbers in this file. 
    # Compare breaker/switch reductions
    for x in reduced_bus_pairs_psse
        _test_psse_reduction_row(x, nr.reverse_bus_search_map)
    end
    #Test number of nonzero elements matches
    @test nnz(Ybus_pnm.data) == length(filter(!iszero, y_values)) + n_ref_bus_elements
    #Test values match 
    for (row_bus, col_bus, val) in zip(row_buses, col_buses, y_values)
        if row_bus ∈ keys(nr.reverse_bus_search_map)
            row_bus = nr.reverse_bus_search_map[row_bus]
        end
        if col_bus ∈ keys(nr.reverse_bus_search_map)
            col_bus = nr.reverse_bus_search_map[col_bus]
        end
        @test isapprox(Ybus_pnm[row_bus, col_bus], val, atol = 1e-3)
    end
end

@testset "14 bus system - re-exported from PowerFlows.jl" begin
    sys = System(joinpath(pwd(), "test", "test_data", "exported_modified_case14_sys 1.raw"))
    Ybus_pnm = Ybus(sys)
    ref_bus_numbers = [
        get_number(x) for
        x in get_components(x -> get_bustype(x) == PSY.ACBusTypes.REF, ACBus, sys)
    ]
    skip_indices = indexin(ref_bus_numbers, Ybus_pnm.axes[1])
    n_ref_bus_elements = nnz(Ybus_pnm.data[skip_indices, :])
    nr = Ybus_pnm.network_reduction

    #Test adjacency and Ybus have same non-zero elements
    findnz(Ybus_pnm.data)[1] == findnz(Ybus_pnm.adjacency_data)[1]
    findnz(Ybus_pnm.data)[2] == findnz(Ybus_pnm.adjacency_data)[2]

    # Compare with PSSE
    Ybus_psse, b_ix_psse, row_buses, col_buses, y_values, reduced_bus_pairs_psse =
        parse_psse_ybus("test/test_data/exported_modified_case14_sys_Ymx.txt")
    # Compare breaker/switch reductions
    for x in reduced_bus_pairs_psse
        _test_psse_reduction_row(x, nr.reverse_bus_search_map)
    end
    #Test number of nonzero elements matches
    @test nnz(Ybus_pnm.data) == length(filter(!iszero, y_values)) + n_ref_bus_elements
    #Test values match 
    for (row_bus, col_bus, val) in zip(row_buses, col_buses, y_values)
        if row_bus ∈ keys(nr.reverse_bus_search_map)
            row_bus = nr.reverse_bus_search_map[row_bus]
        end
        if col_bus ∈ keys(nr.reverse_bus_search_map)
            col_bus = nr.reverse_bus_search_map[col_bus]
        end
        @test isapprox(Ybus_pnm[row_bus, col_bus], val, atol = 1e-3)
    end
end

@testset "WECC 240 bus" begin
    sys_240 = System(
        joinpath(pwd(), "test", "test_data", "240busWECC_2018_PSS33.raw");
        runchecks = false,
    )
    Ybus_pnm = Ybus(sys_240)
    nr = Ybus_pnm.network_reduction
    ref_bus_numbers = [
        get_number(x) for
        x in get_components(x -> get_bustype(x) == PSY.ACBusTypes.REF, ACBus, sys_240)
    ]
    skip_indices = indexin(ref_bus_numbers, Ybus_pnm.axes[1])
    n_ref_bus_elements = nnz(Ybus_pnm.data[skip_indices, :])

    # Compare with PSSE
    Ybus_psse, b_ix_psse, row_buses, col_buses, y_values, reduced_bus_pairs_psse =
        parse_psse_ybus("test/test_data/240busWECC_2018_PSS33_Ymatrix.txt")

    # Compare breaker/switch reductions
    for x in reduced_bus_pairs_psse
        _test_psse_reduction_row(x, nr.reverse_bus_search_map)
    end
    #Test number of nonzero elements matches
    @test nnz(Ybus_pnm.data) == length(filter(!iszero, y_values)) + n_ref_bus_elements
    #Test values match 
    for (row_bus, col_bus, val) in zip(row_buses, col_buses, y_values)
        if row_bus ∈ keys(nr.reverse_bus_search_map)
            row_bus = nr.reverse_bus_search_map[row_bus]
        end
        if col_bus ∈ keys(nr.reverse_bus_search_map)
            col_bus = nr.reverse_bus_search_map[col_bus]
        end
        @test isapprox(Ybus_pnm[row_bus, col_bus], val, atol = 1e-3)
    end
end
