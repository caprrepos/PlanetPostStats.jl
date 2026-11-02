# Time series, multi-year means and standard deviations of 3D dust tracer masses and
# other variables from ROCKE-3D Mars simulations post-processed with scaleacc.
#
# Uses `PlanetPostStats.TimeSeriesMeSd`. All simulation- and user-specific settings are
# the constants in the section "Settings" at the top; no other code needs to be edited.
# Replace the paths, run identifiers, and metadata with your own.

using PlanetPostStats.TimeSeriesMeSd
using OrderedCollections: OrderedDict
using DataFrames

## Settings ===============================================================================

## Simulations ---------------------------------------------------------------------------
const SIMPATH = "/path/to/simulations"      # directory with one subdirectory per run
const POSTDIR = "Post"                      # subdirectory of a run with the scaleacc files
const RUNS = ["MarsDustRun1", "MarsDustRun2"]
const POSTNAMES = ["aij", "aijl", "taijl"]  # scaleacc products read

## Time periods --------------------------------------------------------------------------
const FSIMYEAR = 1        # first year of input data files
const FSIMMONTH = 1       # first month of input data files
const LSIMYEAR = 20       # last year of input data files
const LSIMMONTH = 12      # last month of input data files
const FSTATYEAR = 3       # first year included in time related statistical analysis
const FSTATMONTH = 1      # first month included in time related statistical analysis
const LSTATYEAR = 20      # last year included in time related statistical analysis
const LSTATMONTH = 12     # last month included in time related statistical analysis

## Planet: calendar and seasons ----------------------------------------------------------
const CALENDAR = MARS_CALENDAR
const SEASONS = MARS_SEASONS
const WEIGHT_SEASONS = true    # false: equal weights for the months of a season

## Time coordinate of the output: year TIME_REF_YEAR / month TIME_REF_MONTH gets time = 0,
## i.e. year 1 / month 1 gets time = 1
const TIME_REF_YEAR = 0
const TIME_REF_MONTH = 12
const TIME_UNITS = "months since 0023-12-56"

## Vertical coordinate -------------------------------------------------------------------
const LS1 = 31      # first layer of pressure coordinates above sigma coordinates
const PSF = 6.86    # global mean surface pressure (hPa)
const PTOP = 0.2    # global mean pressure at the bottom of layer LS1 (hPa)

## Optional output on pressure levels: 3D variables also on pressure levels in files
## *.ijp/zmp/gmp.*.nc. `nothing`: no pressure levels; PressureLevels(): default levels of
## the planet (Mars: 43 levels 18.789-0.099 hPa, DEFAULT_PRESSURE_LEVELS["mars"]);
## PressureLevels([7.0, 5.0, 3.0, 1.0]): these levels (hPa); PressureLevels(:model):
## reference levels plm of the input (layer edges Psf, ple); keyword mode=:replace
## writes them instead of the model levels into the regular files.
const PRESSURE_LEVELS = nothing

## Statistics ----------------------------------------------------------------------------
const MISSING_THRESHOLD = 80.0   # percentage threshold required for non-missing data

## Variables -----------------------------------------------------------------------------
const DUST_BINS = ["Clay", "Silt1", "Silt2", "Silt3", "Silt4"]
const DUST_SCALE = 1e-8          # scale of the dust tracers in the scaleacc files
const DUST_UNITS = "kg kg-1"

## Output --------------------------------------------------------------------------------
const OUTPUT_DIR = "Tracermass3d"     # subdirectory of each run for the output
const OUTPUT_NAME = "tracermass3d"    # tag in the output file names

## Metadata ------------------------------------------------------------------------------
const RAW_DATA_ORIGIN = "NASA GISS ROCKE-3D Mars"
const MODEL_REPOSITORY = "<model source repository>"
const GIT_REVISION = "<git revision of the model code>"
const FILE_CREATOR = "<name>"
const ORGANIZATION = "<organization>"
const CONTACT = "<e-mail address>"

## End of settings ========================================================================

## Input variables -----------------------------------------------------------------------
groups = [
    InputGroup("aij", [InputVariable("prsurf"; units="hPa")]),
    InputGroup(
        "aijl",
        [
            InputVariable("z"; units="m"),
            # layer-integrated quantity: stays on model levels with pressure_levels
            InputVariable("airmass"; units="kg m-2", to_pressure=false),
        ],
    ),
    InputGroup(
        "taijl", [InputVariable(b; scale=DUST_SCALE, units=DUST_UNITS) for b in DUST_BINS]
    ),
]

# Accumulated dust mass mixing ratios over the size bins with the names and order of
# PlanetTools.accumulate_dust_masses_modele("", ...): PM2Dust (= Clay), PM4Dust
# (Clay + Silt1), PM8Dust, PM16Dust, PM32Dust (Clay + ... + Silt4).
accumulated = dust_mass_classes(DUST_BINS; units=DUST_UNITS)

## Global attributes of the output files -------------------------------------------------
attributes = OrderedDict{String,Any}(
    "Raw_data_origin" => RAW_DATA_ORIGIN,
    "Model_repository" => MODEL_REPOSITORY,
    "Git_revision" => GIT_REVISION,
    "Input_created_with" => "scaleacc <acc-file> " * join(POSTNAMES, ","),
    "Programming_language" => "Julia v" * string(VERSION),
    "Julia_script_name" => basename(@__FILE__),
    "File_creator" => FILE_CREATOR,
    "Organization" => ORGANIZATION,
    "Contact" => CONTACT,
)

## Configuration -------------------------------------------------------------------------
cfg = TimeSeriesConfig(;
    simpath=SIMPATH,
    postdir=POSTDIR,
    runs=RUNS,
    groups=groups,
    derived=accumulated,
    sim_period=(YearMonth(FSIMYEAR, FSIMMONTH), YearMonth(LSIMYEAR, LSIMMONTH)),
    stat_period=(YearMonth(FSTATYEAR, FSTATMONTH), YearMonth(LSTATYEAR, LSTATMONTH)),
    calendar=CALENDAR,
    seasons=SEASONS,
    time_axis=TimeAxis(;
        reference=YearMonth(TIME_REF_YEAR, TIME_REF_MONTH), units=TIME_UNITS
    ),
    vertical=VerticalCoordinate(LS1, PSF, PTOP),
    pressure_levels=PRESSURE_LEVELS,
    threshold=MISSING_THRESHOLD,
    weight_seasons=WEIGHT_SEASONS,
    output_dir=OUTPUT_DIR,
    output_name=OUTPUT_NAME,
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
