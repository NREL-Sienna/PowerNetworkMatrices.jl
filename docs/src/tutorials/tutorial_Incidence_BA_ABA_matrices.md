# Incidence, BA and ABA matrices

In this tutorial the `IncidenceMatrix`, `BA_matrix` and `ABA_matrix` are presented.
The methods used for their evaluation, as well as how data is stored is shown in
the following subsections.

The matrices here presented are the building blocks for the compuation of the PTDF and LODF matrices.

## IncidenceMatrix

The `PowerNetworkMatrices` package defines the structure `IncidenceMatrix`, which
store the Incidence Matrix of the considered system as well as the most relevant network data.

At first, the `System` data is loaded.

``` @repl tutorial_Incidence_BA_ABA_matrices
using PowerNetworkMatrices
using PowerSystemCaseBuilder

const PNM = PowerNetworkMatrices
const PSB = PowerSystemCaseBuilder

sys = PSB.build_system(PSB.PSITestSystems, "c_sys5");
```

Then the Incidence Matrix is computed as follows:

``` @repl tutorial_Incidence_BA_ABA_matrices
incidence_matrix = PNM.IncidenceMatrix(sys);
```

The `incidence_matrix` variable is a structure of type `IncidenceMatrix`
featuring the following fields:

``` @repl tutorial_Incidence_BA_ABA_matrices
# axis names: row and column names.
# row names: names of the branches
# column names: names of the buses
incidence_matrix.axes

# data: Incidence Matrix
incidence_matrix.data

# lookup: dictionary linking the branche names and bus numbers with the row
# and column numbers, respectively.
incidence_matrix.axes

# ref_bus_positions: set containing the positions of the reference buses.
# this represents the positions where to add the column of zeros. Please refer to the
# exaple in the BA matrix for more details.
incidence_matrix.ref_bus_positions
```

Please note that the matrix data can be easily access by using the following function:

``` @repl tutorial_Incidence_BA_ABA_matrices
PNM.get_data(incidence_matrix)
```

Note that the number of columns is lower than the actual number of system buses since
the column related to the reference bus is discarded.

## BA_Matrix

The `BA_Matrix` is a structure containing the matrix coming from the product of the
`IncidenceMatrix` and the diagonal matrix contianing the impedence of the system's branches ("B" matrix).

The `BA_Matrix` is computed as follows:

``` @repl tutorial_Incidence_BA_ABA_matrices
ba_matrix = PNM.BA_Matrix(sys);
```

As for the `IncidenceMatrix`, here too the `BA_Matrix` structure feature the same fields.

An example related to accessing the matrix data is now provided:
``` @repl tutorial_Incidence_BA_ABA_matrices
# access data by explicitly calling the field "data"
ba_matrix.data

# or by using the "get_data" function
PNM.get_data(ba_matrix)
```
Note that the number of columns is lower than the actual number of system buses since
the column related to the reference bus is discarded.

To add the column of zeros related to the reference bus, it is necessary to use the
information contained in the `ref_bus_positions` field.

``` @repl tutorial_Incidence_BA_ABA_matrices
new_ba_matrix = hcat(
    ba_matrix.data[:,1:collect(ba_matrix.ref_bus_positions)[1]-1],
    zeros(size(ba_matrix, 1), 1),
    ba_matrix.data[:, collect(ba_matrix.ref_bus_positions)[1]:end]
    )
```

However, trying to change the data field with a matrix of different dimension
will result in an error.

``` @repl tutorial_Incidence_BA_ABA_matrices
ba_matrix.data = hcat(
    ba_matrix.data[:,1:collect(ba_matrix.ref_bus_positions)[1]-1],
    zeros(size(ba_matrix, 1), 1),
    ba_matrix.data[:, collect(ba_matrix.ref_bus_positions)[1]:end]
    )
```


## ABA_Matrix

The `ABA_Matrix` is a structure containing the matrix coming from the product of the
`IncidenceMatrix` and the `BA_Matrix`.
It features the same fields as the `IncidenceMatrix` and the `BA_Matrix`, plus the `K` one.
The field `ABA_Matrix.K` stores the LU factorization matrices (using the
methods contained in the package `KLU`).

To evaluate the `ABA_Matrix`, the following command is sufficient:

``` @repl tutorial_Incidence_BA_ABA_matrices
aba_matrix = ABA_Matrix(sys);
```

By default the LU factorization matrices are not computed, leaving the `K` field empty.
In case these are wanted, the keyword `factorize` must be true.

``` @repl tutorial_Incidence_BA_ABA_matrices
aba_matrix_with_LU = ABA_Matrix(sys; factorize=true);

aba_matrix_with_LU.K
```

If the `ABA_Matrix` is already computed but the LU factorization was not performed, this can be done by considering the following command:

``` @repl tutorial_Incidence_BA_ABA_matrices
aba_matrix.K
aba_matrix = factorize(aba_matrix);
aba_matrix.K
```

The following command can then be used to check if the `ABA_Matrix` contains the LU factorization matrices:

``` @repl tutorial_Incidence_BA_ABA_matrices
is_factorized(aba_matrix)
```
