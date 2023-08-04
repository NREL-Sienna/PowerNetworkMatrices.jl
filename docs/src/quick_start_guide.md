# Quick Start Guide

The data for these tutorials is provided in [PowerSystemCaseBuilder](https://github.com/nrel-sienna/PowerSystemCaseBuilder.jl). If you want to build your own case, take a look at the tutorial [Creating and Handling Data for Dynamic Simulations](https://nrel-sienna.github.io/PowerSimulationsDynamics.jl/stable/tutorials/tutorial_dynamic_data/#Creating-and-Handling-Data-for-Dynamic-Simulations).

For more details about loading data and adding more dynamic components check the
[Creating a System with Dynamic devices](https://nrel-sienna.github.io/PowerSystems.jl/stable/modeler_guide/system_dynamic_data/)
section of the documentation in `PowerSystems.jl`.

## Loading data

Data can be loaded from a pss/e raw file and a pss/e dyr file.

``` @repl quick_start_guide
using PowerNetworkMatrices
using PowerSystemCaseBuilder

const PNM = PowerNetworkMatrices
const PSB = PowerSystemCaseBuilder

sys = PSB.build_system(PSB.PSITestSystems, "c_sys5")
```

## Computation of the PTDF matrix

Once system data is loaded, netwrok matrices can be evaluated. The following 
example shows how the PTDF matrix is computed.

The function `PTDF` is called for the evaluation of the matrix and other data. These 
are stored in a structure of type `PTDF`.

``` @repl quick_start_guide
# evaluate the PTDF structure containing the matrix and other data.
ptdf_matrix = PNM.PTDF(sys);

# show the PTDF matrix.
PNM.get_data(ptdf_matrix)
```

As it can be seen, PTDF matrix is stored such that the number of rows is equal 
to the number of buses, number of columns equal to the number of branches.
