using PBWT
using Documenter

DocMeta.setdocmeta!(PBWT, :DocTestSetup, :(using PBWT); recursive=true)

makedocs(;
    modules=[PBWT],
    authors="Roland Laboulaye <rlaboulaye@gmail.com> and contributors",
    repo="https://github.com/rlaboulaye/PBWT.jl/blob/{commit}{path}#{line}",
    sitename="PBWT.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://rlaboulaye.github.io/PBWT.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/rlaboulaye/PBWT.jl",
)
