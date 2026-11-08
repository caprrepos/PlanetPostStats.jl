## Author: Jan P. Perlwitz, Climate, Aerosol, and Pollution Research, LLC
##         (jan.p.perlwitz@caprllc.com)
## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-26 to 2026-11-05

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
  - [`DustOptics`](@ref): optical post-processing of the dust aerosol: dust masses of
    the size bins of the radiative calculations, effective radii, and scattering,
    absorption, and extinction optical depths of the layers and the column at selected
    wavelengths, computed from the `TimeSeriesMeSd` output of the dust tracers and
    processed into the same time series and statistics.

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
include(joinpath("DustOptics", "DustOptics.jl"))

export TimeSeriesMeSd, DustOptics

end # module PlanetPostStats
