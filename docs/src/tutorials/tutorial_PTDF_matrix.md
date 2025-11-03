# PTDF matrix

In this tutorial the methods for computing the Power Transfer Distribution Factors (`PTDF`) are presented.
Before diving into this tutorial we encourage the user to load `PowerNetworkMatrices`, hit the `?` key in the REPL terminal and look for the documention of the different `PTDF` methods available.

## Evaluation of the `PTDF` matrix

The `PTDF` matrix can be evaluated according to two different approaches:
- `Dense`: considers functions for dense matrix multiplication and inversion
- `KLU`: considers functions for sparse matrix multiplication and inversion (**default**)

The evaluation of the `PTDF` matrix can be easily performed starting from importing the system's data and then by simply calling the `PTDF` method.

``` @repl tutorial_PTDF_matrix
using PowerNetworkMatrices
using PowerSystemCaseBuilder

const PNM = PowerNetworkMatrices
const PSB = PowerSystemCaseBuilder

sys = PSB.build_system(PSB.PSITestSystems, "c_sys5");

ptdf_1 = PTDF(sys);

get_ptdf_data(ptdf_1)
```

Advanced users might be interested in computing the `PTDF` matrix starting from either the data contained in the `IncidenceMatrix` and `BA_matrix` structures, or by the information related to the `branches` and `buses` of the system.

``` @repl tutorial_PTDF_matrix
# evaluate the BA_matrix and Incidence_Matrix
ba_matrix = BA_Matrix(sys);
a_matrix = IncidenceMatrix(sys);

# get the PTDF matrix starting from the values of the
# previously computed matrices
ptdf_2 = PTDF(a_matrix, ba_matrix);
get_ptdf_data(ptdf_2)


## Available methods for the computation of the `PTDF` matrix

As previously mentioned, the `PTDF` matrix can be evaluated considering different approaches. The method can be selected by specifying the field `linear_solver` in the `PTDF` function.

``` @repl tutorial_PTDF_matrix
ptdf_dense = PTDF(sys, linear_solver="Dense");
get_ptdf_data(ptdf_dense)

ptdf_klu = PTDF(sys, linear_solver="KLU");
get_ptdf_data(ptdf_klu)
```

By default the "KLU" method is selected, which appeared to require significant less time and memory with respect to "Dense".
Please note that either the `KLU` or `Dense` method is used, the resulting `PTDF` matrix is stored as a dense one.

## Evaluating the `PTDF` matrix considering distributed slack bus

Whenever needed, the `PTDF` matrix can be computed with a distributed slack bus. To do so, a vector of type `Vector{Float64}` needs to be defined, specifying the weights for each bus of the system. These weights identify how the load on the slack bus is redistributed accross the system.

``` @repl tutorial_PTDF_matrix
# consider equal distribution accross each bus for this example
buscount = length(PSY.get_available_components(PSY.ACBus, sys));
dist_slack = 1 / buscount * ones(buscount);
dist_slack_array = dist_slack / sum(dist_slack);
```

Once the vector of the weights is defined, the `PTDF` matrix can be computed by defining the input argument `dist_slack` (empty array `Float64[]` by default):

``` @repl tutorial_PTDF_matrix
ptdf_distr = PTDF(sys, dist_slack=dist_slack_array);
```

The difference between a the matrix computed with and without the `dist_slack` field defined can be seen as follows:

``` @repl tutorial_PTDF_matrix
# with no distributed slack bus
get_ptdf_data(ptdf_klu)
# with distributed slack bus
get_ptdf_data(ptdf_distr)
```

## "Sparse" `PTDF` matrix

The `PTDF` matrix can be computed in a "sparse" fashion by defining the input argument `tol`. If this argument is defined, then elements of the `PTDF` matrix whose absolute values are below the set threshold are dropped. In addition, the matrix will be stored as a sparse one of type `SparseArrays.SparseMatrixCSC{Float64, Int}` instead of `Matrix{Float64}`.

By considering an "extreme" value of 0.2 as `tol`, the `PTDF` matrix can be computed as follows:

``` @repl tutorial_PTDF_matrix
ptdf_sparse = PTDF(sys, tol=0.2);
get_ptdf_data(ptdf_sparse)
```

**NOTE:** 0.2 was used for the purpose of this tutorial. In practice much smaller values are used (e.g., 1e-5).
