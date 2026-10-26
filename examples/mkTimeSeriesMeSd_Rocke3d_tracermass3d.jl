# Time series, multi-year means and standard deviations of 3D dust tracer masses and
# other variables from ROCKE-3D Mars simulations post-processed with scaleacc.
#
# Uses `PlanetPostStats.TimeSeriesMeSd`. All simulation- and user-specific settings are in
# the `TimeSeriesConfig` below; no other code needs to be edited. Replace the paths, run
# identifiers, and metadata with your own.

using PlanetPostStats.TimeSeriesMeSd
using OrderedCollections: OrderedDict
using DataFrames

## Simulations ---------------------------------------------------------------------------
const SIMPATH = "/path/to/simulations"
const RUNS = ["MarsDustRun1", "MarsDustRun2"]
const POSTNAMES = ["aij", "aijl", "taijl"]

## Variables -------------------------------------------------------------------------------
const DUST_BINS = ["Clay", "Silt1", "Silt2", "Silt3", "Silt4"]

groups = [
    InputGroup("aij", [InputVariable("prsurf"; units="hPa")]),
    InputGroup(
        "aijl", [InputVariable("z"; units="m"), InputVariable("airmass"; units="kg m-2")]
    ),
    InputGroup("taijl", [InputVariable(b; scale=1e-8, units="kg kg-1") for b in DUST_BINS]),
]

# Accumulated dust mass mixing ratios over the size bins with the names and order of
# PlanetTools.accumulate_dust_masses_modele("", ...): PM2Dust (= Clay), PM4Dust
# (Clay + Silt1), PM8Dust, PM16Dust, PM32Dust (Clay + ... + Silt4).
accumulated = dust_mass_classes(DUST_BINS; units="kg kg-1")

## Metadata ------------------------------------------------------------------------------
attributes = OrderedDict{String,Any}(
    "Raw_data_origin" => "NASA GISS ROCKE-3D Mars",
    "Model_repository" => "<model source repository>",
    "Git_revision" => "<git revision of the model code>",
    "Input_created_with" => "scaleacc <acc-file> " * join(POSTNAMES, ","),
    "Programming_language" => "Julia v" * string(VERSION),
    "Julia_script_name" => basename(@__FILE__),
    "File_creator" => "<name>",
    "Organization" => "<organization>",
    "Contact" => "<e-mail address>",
)

## Configuration -------------------------------------------------------------------------
cfg = TimeSeriesConfig(;
    simpath=SIMPATH,
    runs=RUNS,
    groups=groups,
    derived=accumulated,
    sim_period=(YearMonth(1, 1), YearMonth(20, 12)),
    stat_period=(YearMonth(3, 1), YearMonth(20, 12)),
    calendar=MARS_CALENDAR,
    seasons=MARS_SEASONS,
    # year 1 / month 1 gets time = 1
    time_axis=TimeAxis(; reference=YearMonth(0, 12), units="months since 0023-12-56"),
    vertical=VerticalCoordinate(31, 6.86, 0.2),
    threshold=80.0,
    weight_seasons=true,       # false: equal weights for the months of a season
    output_dir="Tracermass3d",
    output_name="tracermass3d",
    global_attributes=attributes,
    # filename=(ym, product, run) -> ...,  # only if not "<MON><YYYY>.<product><run>.nc"
)

## Run -----------------------------------------------------------------------------------
results = run_timeseries_statistics(cfg)

## Inspect global means --------------------------------------------------------------------
df = annual_summary(results)
filter(:Varname => ==("prsurf"), df)

df_seasons = seasonal_summary(results)
filter([:Varname, :Season] => (v, s) -> v == "prsurf" && s in ("FMA", "ASO"), df_seasons)

# Optional, with SummaryTables.jl in your environment:
# using SummaryTables
# listingtable(filter(:Varname => ==("prsurf"), df), :MOA; cols=[:Varname], rows=:Experiment)
