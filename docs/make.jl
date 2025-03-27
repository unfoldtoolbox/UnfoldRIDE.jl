using UnfoldRIDE
using Documenter

DocMeta.setdocmeta!(UnfoldRIDE, :DocTestSetup, :(using UnfoldRIDE); recursive = true)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [UnfoldRIDE],
    authors = "Till Prölß, René Skukies, Benedikt Ehinger",
    repo = "https://github.com/unfoldtoolbox/UnfoldRIDE.jl/blob/{commit}{path}#{line}",
    sitename = "UnfoldRIDE.jl",
    format = Documenter.HTML(; canonical = "https://unfoldtoolbox.github.io/UnfoldRIDE.jl"),
    pages = ["index.md"; numbered_pages],
)

deploydocs(; 
    repo = "github.com/unfoldtoolbox/UnfoldRIDE.jl",
    push_preview = true,
    devbranch = "main")
