using SparseArrays

function parse_psse_ybus(path)
    data = readlines(path)
    row_buses, col_buses, y_values = Int[], Int[], ComplexF64[]
    reduced_bus_pairs = []
    values_done = false
    for line in data
        if line == " ZERO IMPEDANCE LINE CONNECTED BUSES (BUS IN MATRIX LISTED FIRST):"
            values_done = true
        end
        if !isempty(line) &&
           line != " ZERO IMPEDANCE LINE CONNECTED BUSES (BUS IN MATRIX LISTED FIRST):"
            if values_done == true
                vals = split(strip(line), r",\s*")
                push!(reduced_bus_pairs, (parse(Int, vals[1]), parse(Int, vals[2])))
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
    unique_buses = sort(unique(vcat(row_buses,col_buses)))
    bus_index = Dict(bus => i for (i, bus) in enumerate(unique_buses))
    row_indices = [bus_index[f] for f in row_buses]
    col_indices = [bus_index[t] for t in col_buses]
    return row_indices, col_indices, y_values, reduced_bus_pairs
end 


@testset "14 bus system" begin 
    sys = System(joinpath(pwd(), "test", "test_data", "modified_14bus_system.raw"))  
    Ybus_pnm = Ybus(sys; check_connectivity=false)
    ref_bus_numbers = [get_number(x) for x in get_components(x->get_bustype(x) == PSY.ACBusTypes.REF, ACBus, sys)]
    skip_indices = indexin(ref_bus_numbers, Ybus_pnm.axes[1])
    n_ref_bus_elements = nnz(Ybus_pnm.data[skip_indices, :])
    nr = Ybus_pnm.network_reduction;

    # Test the NetworkReduction:
    #Test breaker/switch reduction: 
    @test length(nr.bus_reduction_map) == 14
    @test length(nr.reverse_bus_search_map) == 2
    #Test direct branch maps (1:1):
    @test length(keys(nr.direct_branch_map)) == length(keys(nr.reverse_direct_branch_map)) == 13
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
    row_indices, col_indices, y_values, reduced_bus_pairs_psse =
        parse_psse_ybus("test/test_data/modified_14bus_system_Ymx_SWTBRK.txt")  # NOTE: manually changed the ybus numbers in this file. 
    # Compare breaker/switch reductions
    for x in reduced_bus_pairs_psse
        @test (haskey(nr.reverse_bus_search_map, x[1])) ||
            (haskey(nr.reverse_bus_search_map, x[2]))
    end
    #Test number of nonzero elements matches
    @test nnz(Ybus_pnm.data) == length(filter(!iszero, y_values)) + n_ref_bus_elements
    #Test values match 
    for (row_ix, col_ix, val) in zip(row_indices, col_indices, y_values)
        @test isapprox(Ybus_pnm.data[row_ix, col_ix], val, atol = 1e-3)
    end 
end 


@testset "14 bus system - re-exported from PowerFlows.jl" begin 
    sys = System(joinpath(pwd(), "test", "test_data", "exported_modified_case14_sys 1.raw"))  
    Ybus_pnm = Ybus(sys)
    ref_bus_numbers = [get_number(x) for x in get_components(x->get_bustype(x) == PSY.ACBusTypes.REF, ACBus, sys)]
    skip_indices = indexin(ref_bus_numbers, Ybus_pnm.axes[1])
    n_ref_bus_elements = nnz(Ybus_pnm.data[skip_indices, :])
    nr = Ybus_pnm.network_reduction

    #Test adjacency and Ybus have same non-zero elements
    findnz(Ybus_pnm.data)[1] == findnz(Ybus_pnm.adjacency_data)[1]
    findnz(Ybus_pnm.data)[2] == findnz(Ybus_pnm.adjacency_data)[2]

    # Compare with PSSE
    row_indices, col_indices, y_values, reduced_bus_pairs_psse =
        parse_psse_ybus("test/test_data/exported_modified_case14_sys_Ymx.txt")
    # Compare breaker/switch reductions
    for x in reduced_bus_pairs_psse
        @test (haskey(nr.reverse_bus_search_map, x[1])) ||
            (haskey(nr.reverse_bus_search_map, x[2]))
    end
    #Test number of nonzero elements matches
    @test nnz(Ybus_pnm.data) == length(filter(!iszero, y_values)) + n_ref_bus_elements
    #Test values match 
    for (row_ix, col_ix, val) in zip(row_indices, col_indices, y_values)
        @test isapprox(Ybus_pnm.data[row_ix, col_ix], val, atol = 1e-3)
    end 
end 

@testset "WECC 240 bus" begin 
    sys_240 = System(joinpath(pwd(), "test", "test_data", "240busWECC_2018_PSS33.raw"); runchecks= false)
    Ybus_pnm = Ybus(sys_240) 
    ref_bus_numbers = [get_number(x) for x in get_components(x->get_bustype(x) == PSY.ACBusTypes.REF, ACBus, sys_240)]
    skip_indices = indexin(ref_bus_numbers, Ybus_pnm.axes[1])
    n_ref_bus_elements = nnz(Ybus_pnm.data[skip_indices, :])

    # Compare with PSSE
    row_indices, col_indices, y_values, reduced_bus_pairs_psse =
        parse_psse_ybus("test/test_data/240busWECC_2018_PSS33_Ymatrix.txt")

    # Compare breaker/switch reductions
    for x in reduced_bus_pairs_psse
        @test (haskey(nr.reverse_bus_search_map, x[1])) ||
            (haskey(nr.reverse_bus_search_map, x[2]))
    end
    #Test number of nonzero elements matches
    @test nnz(Ybus_pnm.data) == length(filter(!iszero, y_values)) + n_ref_bus_elements
    #Test values match 
    for (row_ix, col_ix, val) in zip(row_indices, col_indices, y_values)
        @test isapprox(Ybus_pnm.data[row_ix, col_ix], val, atol = 1e-3)
    end 
end 
