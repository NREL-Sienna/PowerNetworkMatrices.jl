function brute_force_ptdf(sys::System; tol = 1e-9)
    """
    Calculate the PTDF matrix using a brute-force method.
    """
    data = PowerFlowData(DCPowerFlow(), sys; correct_bustypes = true)
    ptdf = zeros(size(data.branch_activepower_flow_from_to, 1), size(data.bus_magnitude, 1))
    solve_powerflow!(data)
    init_branch_p_from_to = copy(data.branch_activepower_flow_from_to[:, 1])
    for bus_i in 1:size(data.bus_type, 1)
        data.bus_activepower_injection[bus_i, 1] += tol
        solve_powerflow!(data)
        ptdf[:, bus_i] =
            (data.branch_activepower_flow_from_to[:, 1] .- init_branch_p_from_to) ./ tol
        data.bus_activepower_injection[bus_i, 1] -= tol
    end
    return ptdf, data.bus_lookup, data.branch_lookup
end
