# PowerNetworkMatrices.jl

```@meta
CurrentModule = PowerNetworkMatrices
```

## Overview

`PowerNetworkMatrices.jl` is a [`Julia`](http://www.julialang.org) package for
the evaluation of network matrices given the system's data. The package allows to compute
the matrices according to different methods, providing a flexible and powerful tool.

The documentation and code are organized according to the needs of different
users depending on their skillset and requirements. In broad terms there are three categories:

- **Modeler**: Users that want to run a particular analysis or experiment and use `PowerNetworkMatrices.jl` to develop data sets.

- **Model Developer**: Users that want to develop custom components and structs in order to exploit `PowerNetworkMatrices.jl` features to produce custom data sets.

- **Code Base Developers**: Users that want to add new core functionalities or fix bugs in the core capabilities of `PowerNetworkMatrices.jl`.

`PowerNetworkMatrices.jl` is an active project under development, and we welcome your feedback,
suggestions, and bug reports.

## Installation

The latest stable release of PowerNetworkMatrices can be installed using the Julia package manager with

```julia
] add PowerNetworkMatrices
```

For the current development version, "checkout" this package with

```julia
] add PowerNetworkMatrices#main
```

------------
PowerNetworkMatrices has been developed as part of the Scalable Integrated Infrastructure Planning (SIIP) initiative at the U.S. Department of Energy's National Renewable Energy Laboratory ([NREL](https://www.nrel.gov/)).
