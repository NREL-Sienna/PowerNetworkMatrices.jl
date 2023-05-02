sys = System("ACTIVSg2000.m")
sys = System("case_ACTIVSg10k.m")
sys = System("case_ACTIVSg70k.m")

@time begin
    a = IncidenceMatrix(sys)
    ptdf = PTDF(sys)
    lodf = PNM._calculate_LODF_matrix_KLU(a.data, ptdf.data)
end

@time begin
    a = IncidenceMatrix(sys)
    ptdf = PTDF(sys)
    lodf = PNM._calculate_LODF_matrix_DENSE(a.data, ptdf.data)
end

"""
2k
KLU: 0.211567 seconds (71.68 k allocations: 262.962 MiB)
DENSE: 1.929817 seconds (71.64 k allocations: 812.465 MiB, 6.92% gc time)

10k
KLU: 4.986628 seconds (394.14 k allocations: 4.336 GiB, 4.80% gc time, 1.12% compilation time)
DENSE: 104.194184 seconds (316.77 k allocations: 12.754 GiB, 1.25% gc time)

70k

"""
