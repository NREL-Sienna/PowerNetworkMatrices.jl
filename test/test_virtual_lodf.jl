sys = System("ACTIVSg2000.m")
sys = System("case_ACTIVSg10k.m")
sys = System("case_ACTIVSg70k.m"; runchecks = false)

using BenchmarkTools

@benchmark begin
    a = IncidenceMatrix(sys)
    aba = ABA_Matrix(sys; factorize = true)
    ba = BA_Matrix(sys)
    lodf = PNM._calculate_LODF_matrix_KLU(
        a.data,
        aba.K,
        ba.data,
        a.ref_bus_positions,
    )
end

@benchmark begin
    a = IncidenceMatrix(sys)
    ptdf = PTDF(sys)
    lodf_t_1 = PNM._calculate_LODF_matrix_KLU(a.data, ptdf.data)
end

"""
2k
BenchmarkTools.Trial: 48 samples with 1 evaluation.
 Range (min … max):  102.081 ms … 108.081 ms  ┊ GC (min … max): 0.00% … 4.06%
 Time  (median):     103.027 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   104.144 ms ±   1.983 ms  ┊ GC (mean ± σ):  1.15% ± 1.76%

        █▄ ▄                                                     
  ▄▁▁▄▆███▄██▄▁▁█▄▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▆▆▁▆▁▄▄▆▁▄▁▁▁▆▁▄ ▁
  102 ms           Histogram: frequency by time          108 ms <

 Memory estimate: 263.19 MiB, allocs estimate: 77702.

 10k A, K, BA
 BenchmarkTools.Trial: 2 samples with 1 evaluation.
 Range (min … max):  3.983 s …    4.346 s  ┊ GC (min … max): 1.67% … 11.98%
 Time  (median):     4.164 s               ┊ GC (median):    7.05%
 Time  (mean ± σ):   4.164 s ± 256.710 ms  ┊ GC (mean ± σ):  7.05% ±  7.29%

  █                                                        █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  3.98 s         Histogram: frequency by time         4.35 s <

 Memory estimate: 4.34 GiB, allocs estimate: 492284.

 10k PTDF
 BenchmarkTools.Trial: 2 samples with 1 evaluation.
 Range (min … max):  3.354 s …    3.692 s  ┊ GC (min … max): 2.14% … 14.69%
 Time  (median):     3.523 s               ┊ GC (median):    8.72%
 Time  (mean ± σ):   3.523 s ± 238.982 ms  ┊ GC (mean ± σ):  8.72% ±  8.87%

  █                                                        █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  3.35 s         Histogram: frequency by time         3.69 s <

 Memory estimate: 4.33 GiB, allocs estimate: 346900.

"""

@benchmark begin
    a = IncidenceMatrix(sys)
    ptdf = PTDF(sys)
    lodf_t_3 = _calculate_LODF_matrix_KLU2(a.data, ptdf.data)
end

"""
2k
BenchmarkTools.Trial: 36 samples with 1 evaluation.
 Range (min … max):  139.068 ms … 144.730 ms  ┊ GC (min … max): 0.00% … 3.15%
 Time  (median):     139.526 ms               ┊ GC (median):    0.00%
 Time  (mean ± σ):   140.635 ms ±   1.954 ms  ┊ GC (mean ± σ):  0.83% ± 1.33%

   █▃█                                                           
  ▄███▄▄▄▇▇▄▁▁▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▁▁▇▁▁▄▁▁▄▄▄▁▁▄▄▁▁▁▁▄ ▁
  139 ms           Histogram: frequency by time          145 ms <

 Memory estimate: 262.12 MiB, allocs estimate: 77565.

 10k
 BenchmarkTools.Trial: 2 samples with 1 evaluation.
 Range (min … max):  4.271 s …    4.764 s  ┊ GC (min … max): 1.53% … 11.08%
 Time  (median):     4.518 s               ┊ GC (median):    6.57%
 Time  (mean ± σ):   4.518 s ± 348.358 ms  ┊ GC (mean ± σ):  6.57% ±  6.75%

  █                                                        █  
  █▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁█ ▁
  4.27 s         Histogram: frequency by time         4.76 s <

 Memory estimate: 4.33 GiB, allocs estimate: 346757.
"""

@time begin
    a = IncidenceMatrix(sys)
    ptdf = PTDF(sys)
    lodf_t_2 = PNM._calculate_LODF_matrix_DENSE(a.data, ptdf.data)
end

"""
2k
KLU: 0.211567 seconds (71.68 k allocations: 262.962 MiB)
DENSE: 1.929817 seconds (71.64 k allocations: 812.465 MiB, 6.92% gc time)

10k
KLU: 4.986628 seconds (394.14 k allocations: 4.336 GiB, 4.80% gc time, 1.12% compilation time)
DENSE: 104.194184 seconds (316.77 k allocations: 12.754 GiB, 1.25% gc time)

70k
KLU: 

"""
