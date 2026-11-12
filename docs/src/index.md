```@meta
CurrentModule = PlanetPostStats
```

# PlanetPostStats

Statistical post-processing of the monthly output of the GISS post-processing program
`scaleacc` for NASA GISS ModelE and ROCKE-3D simulations
([source code](https://github.com/caprrepos/PlanetPostStats.jl)).

  - **TimeSeriesMeSd**: monthly, seasonal, and annual time series and multi-year means
    and standard deviations of gridded, zonal-mean, and global-mean fields; optional
    output on pressure levels. Start with the [Overview](timeseriesmesd/overview.md).
  - **DustOptics**: effective radii and scattering, absorption, and extinction optical
    depths of the dust aerosol, computed from the `TimeSeriesMeSd` output of the dust
    tracers. Start with the [Overview](dustoptics/overview.md).
  - [Flowcharts](flowcharts.md) of the data flow and the function calls.

```@autodocs
Modules = [PlanetPostStats]
Pages = ["PlanetPostStats.jl"]
```
