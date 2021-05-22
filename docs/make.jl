using HaploSoup
using Documenter

DocMeta.setdocmeta!(HaploSoup, :DocTestSetup, :(using HaploSoup); recursive=true)

makedocs(;
    modules=[HaploSoup],
    authors="Roland Laboulaye <rlaboulaye@gmail.com> and contributors",
    repo="https://github.com/rlaboulaye/HaploSoup.jl/blob/{commit}{path}#{line}",
    sitename="HaploSoup.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://rlaboulaye.github.io/HaploSoup.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/rlaboulaye/HaploSoup.jl",
)
