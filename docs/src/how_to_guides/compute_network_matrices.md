# How to Compute Network Matrices

This guide shows you how to compute various network matrices for your power system.

## Prerequisites

- PowerNetworkMatrices.jl installed
- A power system model loaded (see [Getting Started](@ref))

## Computing PTDF Matrix

To compute the Power Transfer Distribution Factor (PTDF) matrix:

```julia
using PowerNetworkMatrices
const PNM = PowerNetworkMatrices

# Assuming you have a system loaded
ptdf_matrix = PNM.PTDF(sys)

# Access the matrix data
matrix_data = PNM.get_data(ptdf_matrix)
```

## Computing LODF Matrix

To compute the Line Outage Distribution Factor (LODF) matrix:

```julia
lodf_matrix = PNM.LODF(sys)
matrix_data = PNM.get_data(lodf_matrix)
```

## Computing Virtual PTDF Matrix

For systems where you need virtual representation:

```julia
vptdf_matrix = PNM.VirtualPTDF(sys)
matrix_data = PNM.get_data(vptdf_matrix)
```

## Computing Virtual LODF Matrix

Similarly for virtual LODF:

```julia
vlodf_matrix = PNM.VirtualLODF(sys)
matrix_data = PNM.get_data(vlodf_matrix)
```

## Computing Incidence and BA Matrices

For the fundamental network topology matrices:

```julia
# Incidence matrix
incidence_matrix = PNM.IncidenceMatrix(sys)

# BA matrix (Bus-Admittance)
ba_matrix = PNM.BA_Matrix(sys)

# ABA matrix
aba_matrix = PNM.ABA_Matrix(sys)
```

## Working with Pre-computed Matrices

If you have already computed the incidence and BA matrices, you can use them to compute PTDF:

```julia
# Compute base matrices first
ba_matrix = PNM.BA_Matrix(sys)
a_matrix = PNM.IncidenceMatrix(sys)

# Use them to compute PTDF
ptdf_matrix = PNM.PTDF(a_matrix, ba_matrix)
```

## Next Steps

- Learn about [choosing linear solvers](@ref) for optimal performance
- Understand the [theory behind network matrices](@ref) in the Explanation section
