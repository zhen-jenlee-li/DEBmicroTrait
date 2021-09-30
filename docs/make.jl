using DEBmicroTrait
using Documenter

DocMeta.setdocmeta!(DEBmicroTrait, :DocTestSetup, :(using DEBmicroTrait); recursive=true)

makedocs(;
    modules=[DEBmicroTrait],
    authors="giannamars <giannamars@gmail.com> and contributors",
    repo="https://github.com/giannamars/DEBmicroTrait.jl/blob/{commit}{path}#{line}",
    sitename="DEBmicroTrait.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://giannamars.github.io/DEBmicroTrait.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/giannamars/DEBmicroTrait.jl",
)
