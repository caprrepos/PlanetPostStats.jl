# PlanetPostStats.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://caprrepos.github.io/PlanetPostStats.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://caprrepos.github.io/PlanetPostStats.jl/dev/)
[![Build Status](https://github.com/caprrepos/PlanetPostStats.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/caprrepos/PlanetPostStats.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/caprrepos/PlanetPostStats.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/caprrepos/PlanetPostStats.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/JuliaDiff/BlueStyle)

The package generates time-series from monthly-mean two- or three-dimensional output data of `scaleacc`, which converts the accumulation (`acc`) files of NASA GISS ModelE and ROCKE-3D simulations into NetCDF files (e.g. with `aij`, `aijl`, `taijl` in the name, grouping different variables) and calculates zonal and global mean time series from them as well. It also does statistical post-processing of the time series data. Calculation of the dust optics requires the generation of the time series of three-dimensional tracer masses first.

Processing can be done for multiple simulations at once. The package is organized as submodules with their own specific functionalities. Currently, following submodules are available, and more will be added in the near future:

| Submodule        | Purpose |
|:-----------------|:--------|
| `TimeSeriesMeSd` | generates monthly, seasonal, and annual time series and multi-year means and standard deviations of gridded, zonal-mean, and global-mean fields |
| `DustOptics`     | calculates effective radii and scattering, absorption, and extinction optical depths of dust aerosols, and then time series are automatically generated (using `TimeSeriesMeSd`) as output, requires 3D tracer dust mixing ratios generated with `TimeSeriesMeSd` first |

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/caprrepos/PlanetPostStats.jl")
# or, for a local clone
Pkg.develop(path="/path/to/PlanetPostStats.jl")
# once registered in the General registry of Julia:
Pkg.add("PlanetPostStats")
```

Requires Julia ≥ 1.10. Dependencies: PlanetTools ≥ 0.2.2, NCDatasets, DataFrames,
OrderedCollections, CSV, Dates.

## Quick start

Note: The submodules (or their functions) have to be explicitly imported into the workspace to use their functionalities.

```julia
using PlanetPostStats.TimeSeriesMeSd

cfg = TimeSeriesConfig(;
    simpath="/path/to/simulations",
    runs=["MarsDustRun1"],
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
| `pressure_levels` | optional interpolation of 3D variables to pressure levels (`PressureLevels`, see below) |
| `threshold` | minimum percentage of valid data for every averaging step (default 80) |
| `weight_seasons` | weight months by their length in seasonal means (default `true`) |
| `products` | subset of `:mts, :ats, :sts, :moa, :mos, :mom` |
| `representations_3d` | representations written for 3D variables, subset of `:ij, :zm, :gm` (default all); e.g. `[:zm, :gm]` for only zonal and global means |
| `output_dir`, `output_name`, `output_root` | output location and file-name tag |
| `fill_value`, `deflatelevel`, `global_attributes` | NetCDF output settings |
| `filename` | input file-name function (default `scaleacc_filename`, `<MON><YYYY>.<product><run>.nc`) |

## Dust optics

`PlanetPostStats.DustOptics` calculates the dust effective radius of the integrated dust size distribution and scattering, absorption, and extinction optical depths of dust in the size bins and for the integrated dust masses. Its input are the monthly time series of the dust tracer mass mixing ratios (kg kg⁻¹) and of the air mass of the layers (kg m⁻²) written by `TimeSeriesMeSd`
(`<prefix>.<input_name><run>.ijl.mts.nc` in `<simpath>/<run>/<input_dir>/MonthlyTs/`).

```julia
using PlanetPostStats.TimeSeriesMeSd
using PlanetPostStats.DustOptics

bins = [DustBin("Clay1"; source="Clay", radius=0.132e-6, density=2.5e3),
        DustBin("Clay2"; source="Clay", radius=0.23e-6, density=2.5e3),
        DustBin("Clay3"; source="Clay", radius=0.416e-6, density=2.5e3),
        DustBin("Clay4"; source="Clay", radius=0.766e-6, density=2.5e3),
        DustBin("Silt1"; radius=1.386e-6, density=2.65e3)]
sw = SpectralSet("sw", read_spectral_table("aer_sw_dust.nc", "sw_bands.csv"),
                 [0.55e-6, 0.89e-6])
stats = TimeSeriesConfig(; simpath="/path/to/simulations", runs=["MarsDustRun1"],
    sim_period=(YearMonth(1, 1), YearMonth(20, 12)),
    stat_period=(YearMonth(3, 1), YearMonth(20, 12)), calendar=MARS_CALENDAR,
    seasons=MARS_SEASONS, vertical=VerticalCoordinate(31, 6.86, 0.2),
    output_dir="Traceroptics3d", output_name="traceroptics3d")
cfg = DustOpticsConfig(; stats=stats, bins=bins, spectral=[sw],
    fractions=bin_fractions(bins, "Clay" => [0.009, 0.081, 0.234, 0.676]),
    input_dir="Tracermass3d", input_name="tracermass3d")
results = run_dust_optics(cfg)
```

A complete application (ROCKE-3D Mars dust, shortwave and longwave wavelengths, clay
sub-bin fractions selected by run) is in
`examples/mkPostDiagTimeSeriesMeSd_Rocke3d_traceroptics3d.jl`.

## Output

Per run, in `<simpath>/<run>/<output_dir>/` (or `<output_root>/<run>/<output_dir>/`). The
tag is `ij`, `zm`, `gm` (gridded, zonal mean, global mean), or `ijl`, `zml`, `gml` if the
file contains 3D variables on model levels, or `ijp`, `zmp`, `gmp` if it contains 3D
variables only on pressure levels.

With `representations_3d = [:zm, :gm]`, the 3D variables (on model and on pressure
levels) are written only as zonal and global means; the gridded files then contain only
the 2D variables (tag `ij`), and gridded files that would be empty are not created.

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

### References:

  - Hansen, J. E., and L. D. Travis (1974), Light scattering in planetary atmospheres,
    *Space Sci. Rev.*, 16(4), 527–610, doi:10.1007/BF00168069.
  - Heavens, N. G., et al. (2011), The vertical distribution of dust in the Martian
    atmosphere during northern spring and summer: Observations by the Mars Climate
    Sounder and analysis of zonal average vertical dust profiles, *J. Geophys. Res.*,
    116, E04003, doi:10.1029/2010JE003691.
