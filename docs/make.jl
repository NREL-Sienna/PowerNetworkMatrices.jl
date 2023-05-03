using Documenter, PowerNetworkMatrices
import DataStructures: OrderedDict

pages = OrderedDict(
    "Welcome Page" => "index.md",
    "Quick Start Guide" => "quick_start_guide.md",
    "Tutorials" => "tutorials/intro_page.md",
    "Public API Reference" => "api/public.md",
    "Internal API Reference" => "api/internal.md",
)

makedocs(;
    modules = [PowerNetworkMatrices],
    format = Documenter.HTML(;
        mathengine = Documenter.MathJax(),
        prettyurls = haskey(ENV, "GITHUB_ACTIONS")),
    sitename = "PowerNetworkMatrices.jl",
    authors = "Jose Daniel Lara, Alessandro Castelli, Sourabh Dalvi",
    pages = Any[p for p in pages],
    clean = true,
)

deploydocs(;
    repo = "github.com/NREL-Sienna/PowerNetworkMatrices.jl.git",
    target = "build",
    branch = "gh-pages",
    devurl = "dev",
    push_preview = true,
    forcepush = true,
    versions = ["stable" => "v^", "v#.#"],
)
