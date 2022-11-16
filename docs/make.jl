using Documenter, SIIP-PACKAGE


pages = OrderedDict(
        "Welcome Page" => "index.md",
        "Quick Start Guide" => "quick_start_guide.md",
        "Tutorials" =>  "tutorials/intro_page.md",
        "Public API Reference" => "api/public.md",
        "Internal API Reference" => "api/internal.md"
)

pages["Model Library"] = make_model_library(
     categories = [
        Topology,
        StaticInjection,
        Service,
        Branch,
        PSY.DeviceParameter,
    ],
    exceptions = [PSY.DynamicComponent,
                  PSY.ActivePowerControl,
                  PSY.ReactivePowerControl,
                  PSY.DynamicBranch
                  ],
    manual_additions =
        Dict("Service" => ["Reserves" => "model_library/reserves.md"],
        "StaticInjection" => ["Regulation Device" => "model_library/regulation_device.md"],
        "DynamicInjection" => ["Dynamic Inverter" => "model_library/dynamic_inverter.md",
        "Dynamic Generator" => "model_library/dynamic_generator.md",
        ],
        "StaticInjection" => ["HybridSystem" => "model_library/hybrid_device.md"],
        "Branch" => ["Dynamic Lines" => "model_library/dynamic_branch.md"]
        )
)

# postprocess function to insert md
function insert_md(content)
    m = match(r"APPEND_MARKDOWN\(\"(.*)\"\)", content)
    if !isnothing(m)
        md_content = read(m.captures[1], String)
        content = replace(content, r"APPEND_MARKDOWN\(\"(.*)\"\)" => md_content)
    end
    return content
end

# This code performs the automated addition of Literate - Generated Markdowns. The desired
# section name should be the name of the file for instance network_matrices.jl -> Network Matrices
julia_file_filter = x -> occursin(".jl", x)
folders = Dict(
    "Model Library" => filter(julia_file_filter, readdir("docs/src/model_library")),
    "Modeler Guide" => filter(julia_file_filter, readdir("docs/src/modeler_guide")),
    "Model Developer Guide" => filter(julia_file_filter, readdir("docs/src/model_developer_guide")),
    "Code Base Developer Guide" => filter(julia_file_filter, readdir("docs/src/code_base_developer_guide")),
)

for (section, folder) in folders
    for file in folder
        section_folder_name = lowercase(replace(section, " " => "_"))
        outputdir = joinpath(pwd(), "docs", "src", "$section_folder_name")
        inputfile = joinpath("$section_folder_name", "$file")
        infile_path = joinpath(pwd(), "docs", "src", inputfile)
        outputfile = string("generated_", replace("$file", ".jl" => ""))
        execute = occursin("EXECUTE = TRUE", uppercase(readline(infile_path))) ? true : false
        execute && include(infile_path)
        Literate.markdown(infile_path,
                          outputdir;
                          name = outputfile,
                          credit = false,
                          flavor = Literate.DocumenterFlavor(),
                          documenter = true,
                          postprocess = insert_md,
                          execute = execute)
        subsection = titlecase(replace(split(file, ".")[1], "_" => " "))
        push!(pages[section], ("$subsection" =>  joinpath("$section_folder_name", "$(outputfile).md")))
    end
end

makedocs(
    modules = [PowerSystems, InfrastructureSystems],
    format = Documenter.HTML(prettyurls = haskey(ENV, "GITHUB_ACTIONS"),),
    sitename = "PowerSystems.jl",
    authors = "Jose Daniel Lara, Daniel Thom and Clayton Barrows",
    pages = Any[p for p in pages]
)

deploydocs(
    repo = "github.com/NREL-SIIP/PowerSystems.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "master",
    devurl = "dev",
    push_preview=true,
    versions = ["stable" => "v^", "v#.#"],
)
