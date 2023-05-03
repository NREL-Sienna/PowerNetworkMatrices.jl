# PowerNetworkMatrices.jl

[![Main - CI](https://github.com/NREL-Sienna/PowerNetworkMatrices.jl/actions/workflows/main-tests.yml/badge.svg)](https://github.com/NREL-Sienna/PowerNetworkMatrices.jl/actions/workflows/main-tests.yml)
[![codecov](https://codecov.io/gh/NREL-Sienna/PowerNetworkMatrices.jl/branch/main/graph/badge.svg?token=2VvekKsf11)](https://codecov.io/gh/NREL-Sienna/PowerNetworkMatrices.jl)
[![Documentation Build](https://github.com/NREL-Sienna/PowerNetworkMatrices.jl/workflows/Documentation/badge.svg?)](https://nrel-sienna.github.io/PowerNetworkMatrices.jl/stable)
[<img src="https://img.shields.io/badge/slack-@SIIP/PNM-blue.svg?logo=slack">](https://join.slack.com/t/nrel-sienna/shared_invite/zt-glam9vdu-o8A9TwZTZqqNTKHa7q3BpQ)
[![PowerNetworkMatrices.jl Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/PowerNetworkMatrices)](https://pkgs.genieframework.com?packages=PowerNetworkMatrices)

`PowerNetworkMatrices.jl` is able to build classic power systems modeling network matrices such as
[Ybus](https://en.wikipedia.org/wiki/Nodal_admittance_matrix), [PTDF](https://www.powerworld.com/WebHelp/Content/MainDocumentation_HTML/Power_Transfer_Distribution_Factors.htm) and LODF.

## Version Advisory

- PowerNetworkMatrices.jl will work with Julia v1.6+.
- PowerNetworkMatrices.jl exports Matrix methods that were available in PowerSystems.jl version 1.0 and have been implemented as a separate package.

## License

PowerNetworkMatrices is released under a BSD [license](https://github.com/NREL/PowerNetworkMatrices.jl/blob/master/LICENSE).
PowerNetworkMatrices has been developed as part of the Scalable Integrated Infrastructure Planning (SIIP)
initiative at the U.S. Department of Energy's National Renewable Energy Laboratory ([NREL](https://www.nrel.gov/)).
