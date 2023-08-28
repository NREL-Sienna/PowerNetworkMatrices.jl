# LODF matrix

In this tutorial the methods for computing the Line Outage Distribution Factor (`LODF`) are presented.
Before diving into this tutorial we encourage the user to load `PowerNetworkMatrices`, hit the `?` key in the REPL terminal and look for the documentiont of the different `LODF` methods avialable.

## Evaluation of the `LODF` matrix

As for the `PTDF` matrix, the `LODF` one can be evaluated according to two different approaches:
- `Dense`: considers functions for dense matrix multiplication and inversion
- `KLU`: considers functions for sparse matrix multiplication  and inversion(**default**)

The evaluation of the `LODF` matrix can be easily performed starting from importing the system's data and then by simply calling the `LODF` method.

``` @repl tutorial_PTDF_matrix
using PowerNetworkMatrices
using PowerSystemCaseBuilder

const PNM = PowerNetworkMatrices
const PSB = PowerSystemCaseBuilder

# get the System data
sys = PSB.build_system(PSB.PSITestSystems, "c_sys5");

# compute the LODF matrix
lodf_1 = LODF(sys);

lodf_2 = LODF(sys, linear_solver="Dense");

# show matrix
get_lodf_data(lodf_1)
```

Advanced users might be interested in computing the `LODF` matrix starting from either the `branches` and `buses` data (`CASE 1`), the `IncidenceMatrix` and `PTDF` structures (`CASE 2`), or by the information related to `IncidenceMatrix`, `BA_Matrix` and `ABA_Matrix` (`CASE 3`).

``` @repl tutorial_PTDF_matrix
# CASE 1

# get the branches and buses data
branches = PNM.get_ac_branches(sys);
buses = PNM.get_buses(sys);

# compute the LODF matrix from branches and buses data
lodf_3 = LODF(branches, buses);

# CASE 2

# get the Incidence and PTDF matrix
a = IncidenceMatrix(sys);
ptdf = PTDF(sys);

# compute LODF matrix with the two netwrok matrices
lodf_4 = LODF(a, ptdf);

# CASE 3

# get the BA and ABA matrices (ABA matrix must include LU factorization 
# matrices)
ba = BA_Matrix(sys);
aba = ABA_Matrix(sys, factorize = true);

# compute LODF matrix with the three netwrok matrices
lodf_5 = LODF(a, aba, ba);
```

**NOTE:** whenever the method `LODF(sys::System)` is used, the methods previously defined for `CASE 1` and `CASE 2` are executed in sequence. Therefore the method `LODF(a::IncidenceMatrix, ptdf::PTDF)` is the default one when evaluating the `LODF` matrix from the `System` data directly.


## Available methods for the computation of the `LODF` matrix

For those methods that either require the evaluation of the `PTDF` matrix, or that execute this evaluation internally, two different approaches casen be used.

As for the `PTDF` matrix, here too the optional argument `linear_solver` can be specified with either `KLU` (for spars matrix calculation) or `Dense` (for sparse matrix calculation).

``` @repl tutorial_PTDF_matrix
lodf_dense = LODF(sys, linear_solver="Dense");
```

**NOTE (1):** by default the "KLU" method is selected, which appeared to require significant less time and memory with respect to "Dense".
Please note that wether the `KLU` or `Dense` method is used, the resultig `LODF` matrix is stored as a dense one.

**NOTE (2):** for the moment, the method `LODF(a::IncidenceMatrix, aba::ABA_Matrix, ba::BA_Matrix)` will take `KLU` as `linear_solver` option.

## "Sparse" `LODF` matrix

The `LODF` matrix can be computed in a "sparse" fashion by defining the input argument `tol`. If this argument is defined, then elements of the `LODF` matrix whose absolute values are below the set threshold are dropped. In addition, the matrix will be stored as a sparse one of type `SparseArrays.SparseMatrixCSC{Float64, Int64}` type instead of `Matrix{Float64}` one.

By considering an "extreme" value of 0.4 as `tol`, the `LODF` matrix can be computed as follows:

``` @repl tutorial_PTDF_matrix
lodf_sparse = LODF(sys, tol=0.4);
get_lodf_data(lodf_sparse)
```

Please consider that 0.4 was used for the purpose of this tutorial. In practice much smaller values are used (e.g., 1e-5).

**NOTE (1):** elements whose absolute values exceed the `tol` argument are removed from the `LODF` matrix *after* this has been computed.

**NOTE (2):** the `tol` argument does not refer to the "sparsification" tolerance of the `PTDF` matrix that is computed in the `LODF` method.

**NOTE (3):** in case the method `LODF(a::IncidenceMatrix, ptdf::PTDF)` is considerd, an error will be thrown whenever the `tol` argument in the `PTDF` structure used as input is different then 1e-15.