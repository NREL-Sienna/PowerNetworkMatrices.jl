# Load system-specific packages before PNM so relevant PNM extensions
# precompile. The Apple Accelerate backend is now built in to PNM and does
# not need a trigger package; only the MKL/Pardiso extension is gated on a
# `using` call.
@static if (Sys.ARCH === :x86_64 || Sys.ARCH === :i686) && !Sys.isapple()
    import Pkg
    Pkg.add("Pardiso")
    using Pardiso
end

using PowerNetworkMatrices

include("PowerNetworkMatricesTests.jl")
run_tests()
