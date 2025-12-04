# load system-specific packages before PNM: ensures relevant PNM extensions will be loaded.
@static if (Sys.ARCH === :x86_64 || Sys.ARCH === :i686) && !Sys.isapple()
    using MKL
    using Pardiso
else
    using AppleAccelerate
end

using PowerNetworkMatrices

include("PowerNetworkMatricesTests.jl")
run_tests()
