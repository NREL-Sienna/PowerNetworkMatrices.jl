# VirtualLODF

The `VirtualLODF` structure follows the same philosofy as the `VirtualPTDF`: it contains rows of the original `LODF` matrix, evaluated and cached on demand.

Refer to the different arguments of the `VirtualLODF` methods by looking at the "Public API Reference" page.

## How the `VirtualLODF` works

The `VirtualLODF` structure retains many of the similarities of the `VirtualPTDF`. However, its computation is more complex and requires some additional data.

Starting from the system data, the `IncidenceMatrix`, `BA_Matrix` and `ABA_Matrix` (with relative LU factorization matrices) are evaluated. The `ABA_Matrix` and `BA_Matrix` are used for the computation of the diagonal elements of the `PTDF` matrix, and this vector is stored in the `VirtualLODF` structure together with the other structures abovementioned.

Once the `VirtualLODF` is initialized, each row of the matrix can be evaluated separately and on user request. The algorithmic procedure is the following:
1. Define the `VirtualPTDF` structure
2. Call any element of the matrix to define and store the relative row as well as showing the selected element

Regarding point 2, if the row has been stored previosly then the desired element is just loaded from the cache and shown.

The flowchart below shows how the `VirtualLODF` is structured and how it works. Examples will be presented in the following sections.

```@raw html
<img src="../assets/VirtualLODF_scheme.png"/>
```

## Initialize `VirtualLODF` and compute/access row/element

As for the `LODF` matrix, at first the `System` data must be loaded. The "RTS GMLC" systems is considered as example:

``` @repl tutorial_VirtualPTDF_matrix
using PowerNetworkMatrices
using PowerSystemCaseBuilder

const PNM = PowerNetworkMatrices;
const PSB = PowerSystemCaseBuilder;

sys = PSB.build_system(PSB.PSISystems, "RTS_GMLC_DA_sys");
```

At this point the `VirtualLODF` is initialized with the following simple command:

``` @repl tutorial_VirtualPTDF_matrix
v_lodf = VirtualLODF(sys);
```

Now, an element of the matrix can be computed by calling the branch name and bus number:

``` @repl tutorial_VirtualPTDF_matrix
el_C31_2_105 = v_lodf["C31-2", "A3"]
```

This element represent the portion flowing on line "A3" now diverted on line "C31-2" as a consequence of its outage.

Alternatively, the number of the branch and bus (corresponding to the number of the PTDF row and column) can be used. In this case the row and column numbers are mapped by the dictonaries contained in the `lookup` field. 

``` @repl tutorial_VirtualPTDF_matrix
row_number = v_lodf.lookup[1]["C31-2"]
col_number = v_lodf.lookup[2]["A3"]
el_C31_2_105_bis = v_lodf[row_number, col_number]
```

**NOTE**: this example was made for the sake of completeness and considering the actual branch names is reccomended.

As previosly mentioned, in order to evalute a single element of the `VirtualLODF`, the entire row related to the selected branch must be considered. For this reason it is cached for later calls.
This is evident by looking at the following example:

``` @repl tutorial_VirtualPTDF_matrix
sys_2k = PSB.build_system(PSB.PSYTestSystems, "tamu_ACTIVSg2000_sys");

v_lodf_2k = VirtualLODF(sys_2k);

# evaluate PTDF row related to branch "ODESSA 2 0  -1001-ODESSA 3 0  -1064-i_1"
@time v_lodf_2k["ODESSA 2 0  -1001-ODESSA 3 0  -1064-i_1", 
                "BRYAN 1 1   -8159-BRYAN 1 0   -8158-i_1"]

# call same element after the row has been stored
@time v_lodf_2k["ODESSA 2 0  -1001-ODESSA 3 0  -1064-i_1", 
                "BRYAN 1 1   -8159-BRYAN 1 0   -8158-i_1"]
```

## "Sparse" `VirtualPTDF`

Sparsification of each row can be achieved in the same fashion as for the `LODF` matrix, by removing those elements whose absolute values is below a certain threshold.

As for the example show for the `LODF` matrix, here to a very high values of 0.4 is considered for the `tol` field. Again, this value is considered just for the sake of this example.

``` @repl tutorial_VirtualPTDF_matrix
# smaller system for the next examples
sys_2 = PSB.build_system(PSB.PSITestSystems, "c_sys5");

v_lodf_dense = VirtualLODF(sys_2);
v_lodf_sparse = VirtualLODF(sys_2, tol=0.4);
```

Let's now evaluate the same row as before and compare the results:
``` @repl tutorial_VirtualPTDF_matrix
original_row = [v_lodf_dense["1", j] for j in v_lodf_dense.axes[2]]
sparse_row = [v_lodf_sparse["1", j] for j in v_lodf_sparse.axes[2]]
```