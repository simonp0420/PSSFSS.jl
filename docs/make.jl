using PSSFSS
using Documenter

makedocs(;
    clean=false,
    modules=[PSSFSS],
    authors="Peter Simon <psimon0420@gmail.com> and contributors",
    repo="https://github.com/simonp0420/PSSFSS.jl/blob/{commit}{path}#L{line}",
    sitename="PSSFSS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://simonp0420.github.io/PSSFSS.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/simonp0420/PSSFSS.jl",
)
