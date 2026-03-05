using UnfoldRIDE
using Documenter
using Glob
using Literate
using Revise

Revise.revise()


GENERATED = joinpath(@__DIR__, "src", "generated")
SOURCE = joinpath(@__DIR__, "literate")

for subfolder ∈ ["explanations", "HowTo", "tutorials", "reference"]
    local SOURCE_FILES = Glob.glob(subfolder * "/*.jl", SOURCE)
    #config=Dict(:repo_root_path=>"https://github.com/unfoldtoolbox/UnfoldRIDE")
    foreach(fn -> Literate.markdown(fn, GENERATED * "/" * subfolder), SOURCE_FILES)

end


DocMeta.setdocmeta!(UnfoldRIDE, :DocTestSetup, :(using UnfoldRIDE); recursive = true)

makedocs(;
    modules = [UnfoldRIDE],
    authors = "René Skukies, Till Prölss, Benedikt Ehinger",
    #repo="https://github.com/unfoldtoolbox/UnfoldSim.jl/blob/{commit}{path}#{line}",
    repo = Documenter.Remotes.GitHub("unfoldtoolbox", "UnfoldRIDE.jl"),
    sitename = "UnfoldRIDE.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://unfoldtoolbox.github.io/UnfoldRIDE.jl",
        edit_link = "main",
        sidebar_sitename = false,
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Installing Julia & UnfoldRIDE.jl" => "installation.md",
        "Tutorials" => [
            "Quickstart" => "literate/tutorials/11-running_ride.md",
            "Simulate test data" => "literate/tutorials/10-data_simulation.md",
        ],
        "Reference" => [
            "RIDE" => "./generated/reference/overview.md",
        ],
        "HowTo" => [
        ],
        "Developer documentation" => "91-developer.md",
        "API / Docstrings" => "api.md",
    ],
)

deploydocs(;
    repo = "github.com/unfoldtoolbox/UnfoldRIDE.jl",
    devbranch = "main",
    push_preview = true,
)