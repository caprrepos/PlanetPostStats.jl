# PlanetPostStats.jl

Statistical post-processing of `scaleacc` output from NASA GISS ModelE and ROCKE-3D
simulations, for one or many simulations at once.

| Submodule        | Purpose |
|:-----------------|:--------|
| `TimeSeriesMeSd` | monthly, seasonal, and annual time series and multi-year means and standard deviations of gridded, zonal-mean, and global-mean fields |

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/caprrepos/PlanetPostStats.jl")
# or, for a local clone
Pkg.develop(path="/path/to/PlanetPostStats.jl")
```

Requires Julia ≥ 1.13. Dependencies: PlanetTools ≥ 0.2
([caprrepos/PlanetTools.jl](https://github.com/caprrepos/PlanetTools.jl), not registered;
installed from GitHub through the `[sources]` entry of `Project.toml`), NCDatasets,
DataFrames, OrderedCollections, Dates.

## Quick start

```julia
using PlanetPostStats.TimeSeriesMeSd

cfg = TimeSeriesConfig(;
    simpath="/home/jan/Data/Rocke3DMars/Simulations/Planet_20240925",
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
annual_summary(results)        # DataFrame of global-mean MOA/SOA of the 2D variables
```

A complete application (ROCKE-3D Mars dust tracers, 3D fields, accumulated dust size
classes, vertical coordinate, metadata) is in
`examples/mkTimeSeriesMeSd_Rocke3d_tracermass3d.jl`.

## Configuration

All settings are fields of `TimeSeriesConfig`; see its docstring for details.

| Field | Content |
|:------|:--------|
| `simpath`, `runs`, `postdir` | input files in `<simpath>/<run>/<postdir>/` (default `postdir = "Post"`) |
| `groups` | `scaleacc` products and variables to read (`InputGroup`, `InputVariable` with scale and units) |
| `derived` | variables computed after reading (`DerivedVariable`, `cumulative_sums`, `dust_mass_classes`) |
| `sim_period`, `stat_period` | months of the input files; window of the multi-year statistics |
| `calendar` | month lengths and CF calendar name (`MARS_CALENDAR`, `EARTH_NOLEAP_CALENDAR`, `PlanetCalendar`) |
| `seasons` | `SeasonSpec`s (`MARS_SEASONS`, `EARTH_SEASONS`) |
| `time_axis` | reference month and units of the time coordinate (`TimeAxis`) |
| `vertical` | hybrid σ–p coordinate metadata `LS1`, `Psf`, `Ptop` (`VerticalCoordinate`) |
| `threshold` | minimum percentage of valid data for every averaging step (default 80) |
| `weight_seasons` | weight months by their length in seasonal means (default `true`) |
| `products` | subset of `:mts, :ats, :sts, :moa, :mos, :mom` |
| `output_dir`, `output_name`, `output_root` | output location and file-name tag |
| `fill_value`, `deflatelevel`, `global_attributes` | NetCDF output settings |
| `filename` | input file-name function (default `scaleacc_filename`, `<MON><YYYY>.<product><run>.nc`) |

## Output

Per run, in `<simpath>/<run>/<output_dir>/` (or `<output_root>/<run>/<output_dir>/`). The
tag is `ij`, `zm`, `gm` (gridded, zonal mean, global mean), or `ijl`, `zml`, `gml` if the
file contains 3D variables.

| Directory      | Files | Content |
|:---------------|:------|:--------|
| `MonthlyTs`    | `JAN0001-DEC0020.<name><run>.<tag>.mts.nc` | monthly means |
| `AnnualTs`     | `ANN_0001-0020.<name><run>.<tag>.ats.nc` | annual means of complete years |
| `SeasonalTs`   | `FMA_0001-0020.<name><run>.<tag>.sts.nc` | seasonal means, one file per season |
| `AnnualMeSd`   | `ANN_0003_0020.<name><run>.<tag>.{moa,soa}.nc` | multi-year mean / SD of annual means |
| `SeasonalMeSd` | `FMA_0003_0020.<name><run>.<tag>.{mos,sos}.nc` | multi-year mean / SD of seasonal means |
| `MonthlyMeSd`  | `JAN-DEC_0003_0020.<name><run>.<tag>.{mom,som}.nc` | multi-year mean / SD of each month |

Time-series files contain `time` (months from `time_axis.reference` to the first month of
each period) with the CF attributes `calendar` (`"mars"`, `"noleap"`, …) and
`month_lengths`, the integer variables `year(time)` and `month(time)`, and the date
strings `Dates(time)`. `PlanetTools.readnetcdf` returns `time` as numbers. Missing data
are stored as `_FillValue` (`fill_value`, default `-1f30`).

## Statistics

All statistics are computed with `PlanetTools.PlanetStatistics`:

* zonal mean: equal weights over the valid longitudes; global mean: area-weighted;
* annual and seasonal means: month-length weight × availability of the monthly input;
  valid if at least `threshold` % of the months are valid;
* multi-year mean and standard deviation: availability-weighted mean and sample standard
  deviation, `Σw(x − μ)² / (V₁ − V₂/V₁)` with `V₁ = Σw`, `V₂ = Σw²` (`n − 1` for equal
  weights);
* spatial means: valid if at least `threshold` % of the longitudes (zonal) or of the area
  (global) are valid.

One run is held in memory at a time; for each variable all products are computed and
written before the next variable is processed.

## Tests

```julia
using Pkg
Pkg.test("PlanetPostStats")   # synthetic scaleacc files with analytic solutions
```

## Formatting

Code follows Blue Style; format with `using JuliaFormatter; format(".")`
(`.JuliaFormatter.toml` sets `style = "blue"`).
