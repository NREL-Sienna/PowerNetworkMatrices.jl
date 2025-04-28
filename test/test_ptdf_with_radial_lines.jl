@testset "Test PTDF matrix with radial lines" begin
    for name in ["c_sys14", "test_RTS_GMLC_sys"]
        # load the system
        sys = PSB.build_system(PSB.PSITestSystems, name)

        # get the NetworkReduction struct
        rb = get_radial_reduction(sys)

        # get the A and BA matrices without radial lines
        A_rad = IncidenceMatrix(sys; network_reduction = rb)
        BA_rad = BA_Matrix(sys; network_reduction = rb)

        # get the original and reduced PTDF matrices (consider different methods)
        ptdf = PTDF(sys)
        ptdf_rad = PTDF(sys; network_reduction = rb)
        ptdf_rad_A_BA = PTDF(A_rad, BA_rad)

        # check if the same angles and flows are coputed with the matrices of the reduced systems
        # get the indices for the reduced system
        bus_numbers = []
        for i in keys(rb.bus_reduction_map)
            append!(bus_numbers, collect(rb.bus_reduction_map[i]))
        end
        bus_idx = setdiff(
            1:size(ptdf.data, 1),
            append!([ptdf.lookup[1][i] for i in bus_numbers]),
        )
        br_idx =
            setdiff(1:size(ptdf.data, 2), [ptdf.lookup[2][i] for i in rb.removed_branches])

        # now get the injections from the system
        n_buses = length(axes(ptdf, 1))
        bus_lookup = ptdf.lookup[1]
        branch_flow_values = zeros(Float64, length(axes(ptdf, 2)))
        bus_activepower_injection = zeros(Float64, n_buses)
        sources =
            PSY.get_components(d -> !isa(d, PSY.ElectricLoad), PSY.StaticInjection, sys)
        for source in sources
            !PSY.get_available(source) && continue
            bus = PSY.get_bus(source)
            bus_ix = bus_lookup[PSY.get_number(bus)]
            bus_activepower_injection[bus_ix] += PSY.get_active_power(source)
        end
        bus_activepower_withdrawals = zeros(Float64, n_buses)
        loads = PSY.get_components(x -> !isa(x, PSY.FixedAdmittance), PSY.ElectricLoad, sys)
        for l in loads
            !PSY.get_available(l) && continue
            bus = PSY.get_bus(l)
            bus_ix = bus_lookup[PSY.get_number(bus)]
            bus_activepower_withdrawals[bus_ix] += PSY.get_active_power(l)
        end
        power_injection =
            deepcopy(bus_activepower_injection - bus_activepower_withdrawals)

        # get the flows with the PTDF matrix
        ref_flow_values = transpose(ptdf.data) * power_injection

        # evalaute according to the matrix with no radial branches
        reduce_flow_values = zeros((length(br_idx),))
        # change power injection for affrefated leaf buses
        power_injection2 = deepcopy(power_injection)
        for i in keys(rb.bus_reduction_map)
            for j in rb.bus_reduction_map[i]
                power_injection2[ptdf.lookup[1][i]] += power_injection[ptdf.lookup[1][j]]
            end
        end
        reduced_flow_values = transpose(ptdf_rad.data) * power_injection2[bus_idx]

        # now check if flows are the same
        @test isapprox(ref_flow_values[br_idx], reduced_flow_values)
        # for the PTDF from A and BA matrices just need to check the elements
        @test isapprox(ptdf_rad.data, ptdf_rad_A_BA.data, atol = 1e-5)
    end
end

@testset "Test PTDF with radial lines and distributed slack" begin
    for name in ["c_sys14", "test_RTS_GMLC_sys"]
        # load the system
        sys = PSB.build_system(PSB.PSITestSystems, name)
        # get the radial branches
        A = IncidenceMatrix(sys)
        rb = get_radial_reduction(sys)
        # get number of buses
        buscount = length(PNM.get_buses(sys))
        # now compute distributed slack vector
        dist_slack = 1 / buscount * ones(buscount)
        slack_array = dist_slack / sum(dist_slack)
        # adjust to have the same vector with and without leaf node reduction
        bus_numbers = reduce(
            vcat,
            [collect(rb.bus_reduction_map[i]) for i in keys(rb.bus_reduction_map)],
        )
        bus_idx = setdiff(
            1:size(A.data, 2),
            append!([A.lookup[2][i] for i in bus_numbers]),
        )
        br_idx = setdiff(1:size(A.data, 1), [A.lookup[1][i] for i in rb.removed_branches])
        for i in keys(rb.bus_reduction_map)
            for j in rb.bus_reduction_map[i]
                slack_array[A.lookup[2][i]] += slack_array[A.lookup[2][j]]
                slack_array[A.lookup[2][j]] = -9999
            end
        end
        # redefine slack_array
        slack_array[slack_array .== -9999] .= 0
        slack_array1 = slack_array[slack_array .!= -9999]
        # now get the PTDF matrices
        ptdf = PTDF(sys; dist_slack = slack_array)
        ptdf_rad = PTDF(sys; network_reduction = rb, dist_slack = slack_array1)

        # now get the injections from the system
        n_buses = length(axes(ptdf, 1))
        bus_lookup = ptdf.lookup[1]
        branch_flow_values = zeros(Float64, length(axes(ptdf, 2)))
        bus_activepower_injection = zeros(Float64, n_buses)
        sources =
            PSY.get_components(d -> !isa(d, PSY.ElectricLoad), PSY.StaticInjection, sys)
        for source in sources
            !PSY.get_available(source) && continue
            bus = PSY.get_bus(source)
            bus_ix = bus_lookup[PSY.get_number(bus)]
            bus_activepower_injection[bus_ix] += PSY.get_active_power(source)
        end
        bus_activepower_withdrawals = zeros(Float64, n_buses)
        loads = PSY.get_components(x -> !isa(x, PSY.FixedAdmittance), PSY.ElectricLoad, sys)
        for l in loads
            !PSY.get_available(l) && continue
            bus = PSY.get_bus(l)
            bus_ix = bus_lookup[PSY.get_number(bus)]
            bus_activepower_withdrawals[bus_ix] += PSY.get_active_power(l)
        end
        power_injection =
            deepcopy(bus_activepower_injection - bus_activepower_withdrawals)

        # get the flows with the PTDF matrix
        ref_flow_values = transpose(ptdf.data) * power_injection

        # evalaute according to the matrix with no radial branches
        reduce_flow_values = zeros((length(br_idx),))
        # change power injection for affrefated leaf buses
        power_injection2 = deepcopy(power_injection)
        for i in keys(rb.bus_reduction_map)
            for j in rb.bus_reduction_map[i]
                power_injection2[ptdf.lookup[1][i]] += power_injection[ptdf.lookup[1][j]]
            end
        end
        reduced_flow_values = transpose(ptdf_rad.data) * power_injection2[bus_idx]

        # now check if flows are the same
        @test isapprox(ref_flow_values[br_idx], reduced_flow_values)
    end
end

@testset "Test PTDF errors" begin
    # load the system
    sys = PSB.build_system(PSB.PSITestSystems, "c_sys14")

    # get the A and BA matrices without radial lines
    A = IncidenceMatrix(sys)
    BA = BA_Matrix(sys)
    rb = get_radial_reduction(sys)
    BA_rad = BA_Matrix(sys; network_reduction = rb)

    test_value = false
    try
        ptdf = PTDF(A, BA_rad)
    catch err
        if err isa Exception
            test_value = true
        end
    end
    @test test_value
end
