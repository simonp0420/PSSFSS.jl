import Pkg
Pkg.activate(@__DIR__)

using Documenter, DocumenterCitations, DemoCards, PSSFSS

olddir = pwd()
cd(@__DIR__)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:numeric
)

#=
using Literate
cd("literate") do
  include("literate/compile.jl") 
end
=#

demopage, postprocess_cb, demo_assets = makedemos("PSS_&_FSS_Element_Gallery")

assets = String[]
isnothing(demo_assets) || (push!(assets, demo_assets))
push!(assets, "assets/citations.css")

makedocs(;
    clean=false,
    checkdocs = :none,
    warnonly = :cross_references,
    modules=[PSSFSS],
    authors="Peter Simon <psimon0420@gmail.com> and contributors",
    sitename="PSSFSS.jl",
    plugins = [bib,],
    format=Documenter.HTML(;
        repolink="https://github.com/simonp0420/PSSFSS.jl/blob/{commit}{path}#L{line}",
        prettyurls = true,
        canonical="https://simonp0420.github.io/PSSFSS.jl/stable",
        assets=assets,
        size_threshold = nothing,
    ),
    pages=[
        "Contents" => "contents.md",
        "Home" => "index.md",
        "User Manual" => "manual.md",
        demopage,
        "Usage Examples" => "examples.md",
        "API Documentation" => "reference.md",
        "API Index" => "function_index.md",
        "References" => "references.md"
    ],
)
postprocess_cb() # For DemoCards

deploydocs(;
    repo="github.com/simonp0420/PSSFSS.jl",
    devbranch = "main"
)

cd(olddir)