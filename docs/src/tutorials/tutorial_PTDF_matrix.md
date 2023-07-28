# PTDF matrix

In this tutorial the methods for computing the Power Transfer Distribution Factors (`PTDF`) are presented.
Three different methods case be chosen:
- `Dense`: considers functions for dense matrix multiplication and inversion
- `KLU`: considers functions for sparse matrix multiplication and inversion
- `MKLPardiso`: uses MKLPardiso for matrix multiplication and inversion

The abovementioned methods are shortly presented in the following sections.

# Evaluation of the `PTDF` matrix

The evaluation of the PTDF matrix can be easily performed starting from importing the system's data and then by simply calling the `PTDF` function.

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
# previosly cumputed matrices
ptdf_2 = PTDF(a_matrix, ba_matrix);
get_ptdf_data(ptdf_2)

# get the buses and branches of the system
branches = PNM.get_ac_branches(sys);
buses = PNM.get_ac_branches(sys);
ptdf_3 = PTDF(branches, buses);
get_ptdf_data(ptdf_3)
```

NOTE: both the `get_ac_branches` and `get_ac_branches` functions are not exported by the `PowerNetworkMatrices` package, and therefore require the package name to be called as a prefix. However, they are shown here just for the sake of making an example.

# Available methods for the computation of the `PTDF` matrix

As previously mentioned, the `PTDF` matrix can be evaluated considering different approaches. The method can be selected by specifying the field `linear_solver` in the `PTDF` function.

``` @repl tutorial_PTDF_matrix
ptdf_dense = PTDF(sys, linear_solver="Dense")

ptdf_dense = PTDF(sys, linear_solver="MKLPardiso")
```

# Accessing data stored in the `PTDF` structure

# 