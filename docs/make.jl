using StateSpaceRealizables
using Documenter

DocMeta.setdocmeta!(StateSpaceRealizables, :DocTestSetup, :(using StateSpaceRealizables); recursive=true)

makedocs(;
    modules=[StateSpaceRealizables],
    authors="Filip Tronarp <filip.tronarp@uni-tuebingen.de> and contributors",
    repo="https://github.com/filtron/StateSpaceRealizables.jl/blob/{commit}{path}#{line}",
    sitename="StateSpaceRealizables.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://filtron.github.io/StateSpaceRealizables.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/filtron/StateSpaceRealizables.jl",
    devbranch="main",
)
