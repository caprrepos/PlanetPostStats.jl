# PlanetPostStats.jl

Statistical post-processing of `scaleacc` output from NASA GISS ModelE and ROCKE-3D
simulations, for one or many simulations at once.

| Submodule        | Replaces                                                    | Purpose |
|:-----------------|:------------------------------------------------------------|:--------|
| `TimeSeriesMeSd` | `mkTimeSeriesMeSd_*.jl` + `griddedPlanetFields*_{calculate,write}MeSd.jl` | monthly / seasonal / annual time series and multi-year means & standard deviations of gridded, zonal-mean and global-mean fields |

Further scripts are added as further submodules (`src/<Name>/<Name>.jl`, included in
`src/PlanetPostStats.jl`).

## Installation (private repository)

```julia
using Pkg
Pkg.develop(path="/path/to/PlanetPostStats")          # local clone
# or
Pkg.add(url="git@github.com:<user>/PlanetPostStats.jl.git")
```

Requires Julia ≥ 1.13 and PlanetTools ≥ 0.2. Dependencies: PlanetTools (private,
`caprrepos/PlanetTools.jl`, resolved through the `[sources]` entry of `Project.toml`, so
your GitHub credentials must give `git` access to it), NCDatasets, DataFrames,
OrderedCollections, Dates.

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

The complete conversion of `mkTimeSeriesMeSd_Rocke3d_tracermass3d.jl` is in
`examples/mkTimeSeriesMeSd_Rocke3d_tracermass3d.jl`.

## Mapping of the old script constants

| Old constant(s)                                   | `TimeSeriesConfig` field |
|:--------------------------------------------------|:-------------------------|
| `BASEPATH`, `SBRANCH` → `SIMPATH`                 | `simpath` |
| `RUN_ID`                                          | `runs` |
| `POSTDIR`                                         | `postdir` (default `"Post"`) |
| `POSTNAMES`, `PHYSNAMES`, `PHYSSCALES`, `PHYSUNITS` | `groups` (`InputGroup`, `InputVariable`) |
| `accumulate_dust_masses_modelE`                   | `derived` (`dust_mass_classes`, `cumulative_sums`, `DerivedVariable`) |
| `GROUPDIR`, `GROUPNAME`                           | `output_dir`, `output_name` |
| `FSIMYEAR` … `LSIMMONTH`                          | `sim_period` |
| `FSTATYEAR` … `LSTATMONTH`                        | `stat_period` |
| `PLANET`, `DAYS_IN_MONTH`, `DAY_LENGTH`, `CMONLIST` | `calendar` (`MARS_CALENDAR`, `EARTH_NOLEAP_CALENDAR`, `PlanetCalendar`) |
| `FIRST_MONTH_OF_SEASON`, `LENGTH_OF_SEASON`, `SEASON_NAMES` | `seasons` (`SeasonSpec`, `MARS_SEASONS`, `EARTH_SEASONS`) |
| `TIME_UNITS`                                      | `time_axis` (`TimeAxis`) |
| `LS1`, `PSF`, `PTOP`                              | `vertical` (`VerticalCoordinate`) |
| `MISSING_THRESHOLD`                               | `threshold` |
| `FILL_VALUE`                                      | `fill_value` |
| `generalattributes`                               | `global_attributes` |
| `fnames_modelE_2d`                                | `filename` (default `scaleacc_filename` = `PlanetTools.fnames_modele`) |

## Output

Per run, in `<simpath>/<run>/<output_dir>/` (file names as before; tag `ij/zm/gm`, or
`ijl/zml/gml` if the file contains 3D variables):

| Directory      | Files |
|:---------------|:------|
| `MonthlyTs`    | `JAN0001-DEC0020.<name><run>.<tag>.mts.nc` |
| `AnnualTs`     | `ANN_0001-0020.<name><run>.<tag>.ats.nc` |
| `SeasonalTs`   | `FMA_0001-0020.<name><run>.<tag>.sts.nc` (one per season) |
| `AnnualMeSd`   | `ANN_0003_0020.<name><run>.<tag>.{moa,soa}.nc` |
| `SeasonalMeSd` | `FMA_0003_0020.<name><run>.<tag>.{mos,sos}.nc` |
| `MonthlyMeSd`  | `JAN-DEC_0003_0020.<name><run>.<tag>.{mom,som}.nc` |

Time-series files contain `time` (months from `time_axis.reference` to the first month of
each period) with the CF attributes `calendar` (`"mars"`, `"noleap"`, …) and
`month_lengths`, the integer variables `year(time)` and `month(time)`, and the strings
`Dates(time)`. `PlanetTools.readnetcdf` returns `time` as numbers.

## Statistics: PlanetTools conventions

All statistics are computed by `PlanetTools.PlanetStatistics` (`calc_planet_means!`,
`calc_annual_mean_ts!`, `calc_seasonal_mean_ts!`, `calc_multiyear_stats!`) through thin
adapters (`planet_means`, `annual_means`, `seasonal_means`, `multiyear_stats`,
`monthly_stats`) that reshape one variable of one run to PlanetTools' `(Sim, Var, …)`
layout without copying. Hence:

* zonal mean: equal weights over valid longitudes; global mean: area-weighted;
* annual/seasonal means: month-length weight × availability of the monthly input;
  threshold on the *number* of valid months;
* multi-year mean and SD: availability-weighted mean and **sample** SD
  (`Σw(x − μ)² / (V₁ − V₂/V₁)`, i.e. `n − 1` for equal weights);
* one `threshold` (percent) for every step.

## Design decisions

* **One run in memory at a time**; per variable, all products are computed and written
  before the next variable, so peak memory ≈ one run's input + one variable's products.
* **Statistics from PlanetTools**, file names from `fnames_modele`, dust classes from the
  bin/class tables of `accumulate_dust_masses_modele` — one implementation to maintain.
  With the loop kernels of PlanetTools ≥ 0.2 the whole statistics chain of one
  72 × 46 × 40 × 240 variable takes ≈ 0.5 s and 0.13 GiB (2 threads).
* **Output with `PlanetTools.ReadWritePlanet.writenetcdf`** (≥ 0.2), which writes missing
  values as `_FillValue` and the planetary `time` axis without CF date encoding; the time
  coordinate carries the CF attributes `calendar` (e.g. `"mars"`) and `month_lengths`.

## Verification

```julia
include("test/runtests.jl")   # synthetic scaleacc files with analytic solutions
```

## Formatting

Code follows Blue Style; format with `using JuliaFormatter; format(".")`
(`.JuliaFormatter.toml` sets `style = "blue"`).
