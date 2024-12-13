using PSSFSS
using Documenter
using DemoCards

#=
using Literate
cd("literate") do
  include("literate/compile.jl") 
end
=#

demopage, postprocess_cb, demo_assets = makedemos("PSS_&_FSS_Element_Gallery")
@show demopage
assets = String[]
isnothing(demo_assets) || (push!(assets, demo_assets))

makedocs(;
    clean=false,
    checkdocs = :none,
    warnonly = :cross_references,
    modules=[PSSFSS],
    authors="Peter Simon <psimon0420@gmail.com> and contributors",
    sitename="PSSFSS.jl",
    format=Documenter.HTML(;
        repolink="https://github.com/simonp0420/PSSFSS.jl/blob/{commit}{path}#L{line}",
        prettyurls = true,
        canonical="https://simonp0420.github.io/PSSFSS.jl/stable",
        assets=assets,
        size_threshold = nothing,
    ),
    pages=[
        "Home" => "index.md",
        "User Manual" => "manual.md",
        demopage,
        "Usage Examples" => "examples.md",
        "Function Reference" => "reference.md",
        "Index" => "function_index.md"
    ],
)
postprocess_cb() # For DemoCards

deploydocs(;
    repo="github.com/simonp0420/PSSFSS.jl",
    devbranch = "main"
)
