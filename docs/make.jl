using PlanetPostStats
using Documenter
using DocumenterMermaid

DocMeta.setdocmeta!(
    PlanetPostStats, :DocTestSetup, :(using PlanetPostStats); recursive=true
)

makedocs(;
    modules=[PlanetPostStats],
    authors="Jan P Perlwitz <jan.p.perlwitz@caprllc.com>",
    sitename="PlanetPostStats.jl",
    format=Documenter.HTML(;
        canonical="https://caprrepos.github.io/PlanetPostStats.jl",
        edit_link="main",
        assets=String[],
        collapselevel=1,
    ),
    pages=[
        "Home" => "index.md",
        "TimeSeriesMeSd" => [
            "Overview" => "timeseriesmesd/overview.md",
            "Configuration" => "timeseriesmesd/configuration.md",
            "Calendar and seasons" => "timeseriesmesd/calendar.md",
            "Input: scaleacc files" => "timeseriesmesd/input.md",
            "Statistics" => "timeseriesmesd/statistics.md",
            "Pressure levels" => "timeseriesmesd/pressure_levels.md",
            "Output files" => "timeseriesmesd/output.md",
            "Driver" => "timeseriesmesd/driver.md",
        ],
        "DustOptics" => [
            "Overview" => "dustoptics/overview.md",
            "Size bins and kernels" => "dustoptics/bins.md",
            "Optical lookup tables" => "dustoptics/spectral.md",
            "Configuration" => "dustoptics/configuration.md",
            "Input" => "dustoptics/input.md",
            "Driver" => "dustoptics/driver.md",
        ],
        "Flowcharts" => "flowcharts.md",
        "Index" => "api_index.md",
    ],
)

deploydocs(; repo="github.com/caprrepos/PlanetPostStats.jl", devbranch="main")
