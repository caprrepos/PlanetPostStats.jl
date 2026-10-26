"""
    PlanetPostStats

Statistical post-processing of output from the GISS post-processing program `scaleacc`,
which converts the accumulation (`acc`) files of NASA GISS ModelE and ROCKE-3D
simulations into monthly NetCDF files (e.g. the `aij`, `aijl`, `taijl` products).

The package is organized into submodules, one per post-processing task. Each submodule
receives *all* simulation- and user-specific information (paths, run identifiers,
variables, calendars, seasons, time windows, metadata) as function arguments, so that no
source code has to be edited when the post-processing is applied to a different set of
simulations or by a different user.

# Submodules

  - [`TimeSeriesMeSd`](@ref): extraction of monthly, seasonal, and annual time series of
    gridded (`ij`), zonal-mean (`zm`), and global-mean (`gm`) fields from a set of
    simulations, and calculation of multi-year means and standard deviations of annual
    means, seasonal means, and the individual months of the annual cycle within a
    selectable time window. All results are written to NetCDF files.

Further submodules for other post-processing tasks are meant to be added as `src/<SubmoduleName>/<SubmoduleName>.jl` and included below.

# Example

```julia
using PlanetPostStats
using PlanetPostStats.TimeSeriesMeSd

cfg = TimeSeriesConfig(;
    simpath="/path/to/simulations",
    runs=["MarsDustRun1"],
    groups=[InputGroup("aij", [InputVariable("prsurf"; units="hPa")])],
    sim_period=(YearMonth(1, 1), YearMonth(20, 12)),
    stat_period=(YearMonth(3, 1), YearMonth(20, 12)),
    calendar=MARS_CALENDAR,
    seasons=MARS_SEASONS,
    output_dir="Prsurf",
    output_name="prsurf",
)
results = run_timeseries_statistics(cfg)
```
"""
module PlanetPostStats

include(joinpath("TimeSeriesMeSd", "TimeSeriesMeSd.jl"))

export TimeSeriesMeSd

end # module PlanetPostStats
