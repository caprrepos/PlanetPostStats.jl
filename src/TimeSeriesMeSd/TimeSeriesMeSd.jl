"""
    TimeSeriesMeSd

Submodule of `PlanetPostStats` that turns monthly `scaleacc` output of one or more
ModelE/ROCKE-3D simulations into

 1. monthly mean time series (`mts`),
 2. annual mean time series (`ats`) of full calendar years,
 3. seasonal mean time series (`sts`), one series per season,
 4. multi-year means and standard deviations of the annual means (`moa`, `soa`),
 5. multi-year means and standard deviations of the seasonal means (`mos`, `sos`),
 6. multi-year means and standard deviations of every month of the annual cycle
    (`mom`, `som`),

each of them for three spatial representations:

  - `ij`  — the gridded field (longitude × latitude [× level]),
  - `zm`  — the area-weighted zonal mean (latitude [× level]),
  - `gm`  — the area-weighted global mean ([level]).

Two-dimensional (`lon × lat`) and three-dimensional (`lon × lat × level`) variables can be
mixed; they are written into the same files. Files containing at least one 3D variable
get the representation tags `ijl`, `zml`, `gml` (e.g. `*.ijl.mts.nc`), files with only 2D
variables get `ij`, `zm`, `gm`.

# Workflow

All settings that may differ between simulations or users are collected in a
[`TimeSeriesConfig`](@ref). The single entry point [`run_timeseries_statistics`](@ref)
then processes every run listed in the configuration *one run at a time* (bounded
memory), and returns a vector of [`RunResult`](@ref)s containing the paths of the written
files and small summary tables of global means for interactive inspection.

# Treatment of missing data

All statistics are computed with `PlanetTools.PlanetStatistics` (see the adapters in
`statistics.jl`), so the definitions are those of PlanetTools. Missing data are `NaN32`.
One validity threshold `threshold` (percent, `TimeSeriesConfig.threshold`, default 80)
applies to every step:

  - zonal means: equal-weight mean over the valid longitudes, valid if at least
    `threshold` % of the longitudes are valid (availability = valid fraction);
  - global means: area-weighted, valid if at least `threshold` % of the area is valid
    (availability = valid area fraction);
  - annual and seasonal means: month weights (month length, see `weight_seasons`) times
    the availability of the input month; valid if at least `threshold` % of the months of
    the period are valid (availability = fraction of valid months);
  - multi-year statistics (`moa/soa`, `mos/sos`, `mom/som`): availability-weighted mean
    and *sample* standard deviation (`Σw(x − μ)² / (V₁ − V₂/V₁)`, i.e. `n − 1` for equal
    weights), weights = availabilities of the annual/seasonal means (or of the monthly
    zonal/global means); valid if at least `threshold` % of the samples are valid; the
    standard deviation needs at least two valid samples.

Output files are written with `PlanetTools.ReadWritePlanet.writenetcdf`. Invalid results
are `NaN32` in memory and the fill value (`TimeSeriesConfig.fill_value`,
default `-1f30`) in the output files, with a matching `_FillValue` attribute.

# Example

```julia
using PlanetPostStats.TimeSeriesMeSd

cfg = TimeSeriesConfig(;
    simpath="/data/Rocke3DMars/Simulations/Planet_20240925",
    runs=["PS20240925TdustM40_Mars_AMF0p5w"],
    groups=[
        InputGroup("aij", [InputVariable("prsurf"; units="hPa")]),
        InputGroup("taijl", [InputVariable("Clay"; scale=1e-8, units="kg kg-1")]),
    ],
    sim_period=(YearMonth(1, 1), YearMonth(20, 12)),
    stat_period=(YearMonth(3, 1), YearMonth(20, 12)),
    calendar=MARS_CALENDAR,
    seasons=MARS_SEASONS,
    output_dir="Tracermass3d",
    output_name="tracermass3d",
)
results = run_timeseries_statistics(cfg)
annual_summary(results)
```
"""
module TimeSeriesMeSd

using Dates: Dates, Date, UTC, Year
using NCDatasets: NCDatasets, NCDataset, dimnames
using OrderedCollections: OrderedDict
using DataFrames: DataFrame
using PlanetTools.PlanetStatistics:
    calc_planet_means!, calc_annual_mean_ts!, calc_seasonal_mean_ts!, calc_multiyear_stats!
using PlanetTools.Constructors: fnames_modele
using PlanetTools.ReadWritePlanet: writenetcdf
using PlanetTools.DerivePhysicsVariablesModelE: DerivePhysicsVariablesModelE

include("calendar.jl")
include("config.jl")
include("statistics.jl")
include("scaleacc_input.jl")
include("netcdf_output.jl")
include("driver.jl")

# Calendar and time handling
export YearMonth, PlanetCalendar, SeasonSpec, TimeAxis
export MARS_CALENDAR, EARTH_NOLEAP_CALENDAR, MARS_SEASONS, EARTH_SEASONS
export yearmonth_range, months_between, month_weights, days_in_year, monthname
export Period, season_months, time_value
export annual_periods, seasonal_periods, climatology_slots, period_in_window

# Configuration
export InputVariable, InputGroup, DerivedVariable, VerticalCoordinate, TimeSeriesConfig
export ALL_PRODUCTS, sum_of, cumulative_sums, dust_mass_classes, scaleacc_filename
export validate_config
export output_directory

# Numerical kernels
# Adapters to PlanetTools.PlanetStatistics
export SpatialMeans, PeriodMeans
export planet_means, annual_means, seasonal_means, multiyear_stats, monthly_stats

# I/O and driver
export load_run, RunData, GridInfo, VariableData
export run_timeseries_statistics, process_run, RunResult, annual_summary, seasonal_summary

end # module TimeSeriesMeSd
