using PlanetPostStats
using Documenter

DocMeta.setdocmeta!(PlanetPostStats, :DocTestSetup, :(using PlanetPostStats); recursive=true)

makedocs(;
    modules=[PlanetPostStats],
    authors="Jan P Perlwitz <jan.p.perlwitz@caprllc.com>",
    sitename="PlanetPostStats.jl",
    format=Documenter.HTML(;
        canonical="https://caprrepos.github.io/PlanetPostStats.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/caprrepos/PlanetPostStats.jl",
    devbranch="main",
)
