
# Extensions are loaded when trigger packages (Pardiso, AppleAccelerate) are loaded

# Check if MKL/Pardiso extension is available at runtime
function _has_mkl_pardiso_ext()
    ext = Base.get_extension(@__MODULE__, :MKLPardisoExt)
    return !isnothing(ext)
end

_mkl_pardiso_install_error() =
    """The MKL/Pardiso extension is not available.
    Install the Pardiso package:
    julia> using Pkg; Pkg.add(\"Pardiso\")"""

# Note that MKL can be used on Apple macOS x86_64, however the user will
# need to pin MKL_jll to version 2023 in their project environment:
# using Pkg
# if Sys.isapple() && (Sys.ARCH == :x64_64)
#     Pkg.add(name="MKL_jll"; version = "2023")
# end

# Check if AppleAccelerate extension is available at runtime
function _has_apple_accelerate_ext()
    ext = Base.get_extension(@__MODULE__, :AppleAccelerateExt)
    return !isnothing(ext)
end

_apple_accelerate_install_error() =
    """The AppleAccelerate extension is not available.
    This solver is only available on macOS.
    Install AppleAccelerate:
    julia> using Pkg; Pkg.add(\"AppleAccelerate\")"""

# _create_apple_accelerate_factorization is defined in ext/AppleAccelerateExt.jl
# when AppleAccelerate package is loaded

"Set a preference of the backend library for linear algebra operations."
function set_linalg_backend_preference(linalglib::Union{String, Nothing})
    if !isnothing(linalglib) && !(linalglib in SUPPORTED_LINEAR_SOLVERS)
        throw(ArgumentError("Unsupported linear algebra backend requested: $(linalglib)"))
    end
    Preferences.@set_preferences!("linalg_backend" => linalglib)
    @info("""Linear algebra backend library: preference set to $(linalglib);
          you may need to restart your Julia session for this change to take effect.""")
end

set_linalg_backend_preference() = set_linalg_backend_preference(nothing)
set_linalg_backend_preference(linalglib::Symbol) =
    set_linalg_backend_preference(String(linalglib))

get_linalg_backend_preference() = Preferences.@load_preference("linalg_backend")

"Set a preference whether to run check_linalg_backend at the package loading time."
set_linalg_backend_check(check::Bool) =
    Preferences.@set_preferences!("linalg_backend_check" => check)

get_linalg_backend_check() = Preferences.@load_preference("linalg_backend_check")

"""
Check for the recommended linear algebra backend for the
user's operating system and whether it matches the preferred one.
"""
function check_linalg_backend()
    user_linalg_backend = get_linalg_backend_preference()
    if !isnothing(user_linalg_backend)
        @info """The linear algebra library preference has been set to $(user_linalg_backend).
                To change this for your active project, call the function
                PowerNetworkMatrices.set_linalg_backend_preference()
                with one of "MKLPardiso" or "AppleAccelerate", or `nothing` to turn off.
              """
    end

    no_msg(lib) = """For faster dense matrix operations, consider using $(lib):
                  pkg> add $(lib) # if not in your active project
                  using $(lib)    # after loading PowerNetworkMatrices and before any matrix operations
                Sparse factorization still uses KLU (recommended)."""
    go_msg(lib) = "The linear algebra backend $(lib) is loaded."
    yo_msg(lib) = """The $(lib) extension for PowerNetworkMatrices is not loaded
                     even though the corresponding linear algebra library was requested."""
    if !Sys.isapple() || user_linalg_backend == "MKLPardiso"
        if _has_mkl_pardiso_ext()
            @info go_msg("MKLPardiso")
        else
            if user_linalg_backend == "MKLPardiso"
                @warn yo_msg("MKLPardiso")
            end
            @info no_msg("Pardiso")
            @info "See https://github.com/JuliaSparse/Pardiso.jl for more details."
        end
    elseif Sys.isapple()
        if _has_apple_accelerate_ext()
            @info go_msg("AppleAccelerate")
        else
            if user_linalg_backend == "AppleAccelerate"
                @warn yo_msg("AppleAccelerate")
            end
            @info no_msg("AppleAccelerate")
            @info "See https://github.com/JuliaLinearAlgebra/AppleAccelerate.jl"
        end
    end
end
