# Time series, multi-year means and standard deviations of 3D dust tracer masses and
# other variables from ROCKE-3D Mars simulations post-processed with scaleacc.
#
# Uses `PlanetPostStats.TimeSeriesMeSd`. All simulation- and user-specific settings are
# the constants in the section "Settings" at the top; no other code needs to be edited.
# Replace the paths, run identifiers, and metadata with your own.

using DataFrames
using OrderedCollections: OrderedDict
using PlanetPostStats.TimeSeriesMeSd

## Settings ===============================================================================

## Simulations ---------------------------------------------------------------------------
const SIMPATH = "/path/to/simulations"      # directory with one subdirectory per run
const POSTDIR = "Post"                      # subdirectory of a run with the scaleacc files
const RUNS = ["Run1", "Run2"]
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

## Variables -----------------------------------------------------------------------------
const DUST_BINS = ["Clay", "Silt1", "Silt2", "Silt3", "Silt4"]
const DUST_SCALE = 1e-8          # scale of the dust tracers in the scaleacc files
const DUST_UNITS = "kg kg-1"

## Statistics ----------------------------------------------------------------------------
const MISSING_THRESHOLD = 80.0   # percentage threshold required for non-missing data
const FILL_VALUE = -1.0f30

## Output --------------------------------------------------------------------------------
const OUTPUT_DIR = "Tracermass3d"     # subdirectory of each run for the output
const OUTPUT_NAME = "tracermass3d"    # tag in the output file names

## Metadata ------------------------------------------------------------------------------
const RAW_DATA_ORIGIN = "NASA GISS ROCKE-3D Mars"
const MODEL_REPOSITORY = "<model source repository>"
const MODEL_GIT_REVISION = "<git revision of the model code>"
const FILE_CREATOR = "<name>"
const ORGANIZATION = "<organization>"
const CONTACT = "<e-mail address>"

## Global means shown after the run ----------------------------------------------------
const INSPECT_VARIABLES = ["prsurf"]
const INSPECT_SEASONS = ["FMA", "MJJ", "ASO", "NDJ"]

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
    "Model_git_revision" => MODEL_GIT_REVISION,
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
    fill_value=FILL_VALUE,
    output_dir=OUTPUT_DIR,
    output_name=OUTPUT_NAME,
    global_attributes=attributes,
    # filename=(ym, product, run) -> ...,  # only if not "<MON><YYYY>.<product><run>.nc"
)

## Read, calculate, and write ------------------------------------------------------------
# Per run: monthly, annual, and seasonal time series and multi-year means and standard
# deviations (annual, seasonal, monthly) of gridded, zonal-mean, and global-mean fields,
# written to <SIMPATH>/<run>/<OUTPUT_DIR>/{MonthlyTs,AnnualTs,SeasonalTs,AnnualMeSd,
# SeasonalMeSd,MonthlyMeSd}.
results = run_timeseries_statistics(cfg)
[(r.run, length(r.months), length(r.files)) for r in results]

## Inspect global means --------------------------------------------------------------------
# Long-time mean of annual global means
df = annual_summary(results)
# Long-time mean of seasonal global means
df_seasons = seasonal_summary(results)

# Show as tables
# Optional, with SummaryTables.jl in your environment:
#using SummaryTables
#println("Long-time mean of annual global means")
#for var in INSPECT_VARIABLES
#    dfsel = filter(:Varname => v -> v == var, df)
#    display(listingtable(dfsel, :MOA; cols=[:Varname], rows=:Experiment, sort=false))
#end

#println("Long-time mean of seasonal global means")
#for var in INSPECT_VARIABLES
#    dfsel = filter([:Varname, :Season] => (v, s) -> v == var && s in INSPECT_SEASONS,
#        df_seasons)
#    display(listingtable(dfsel, :MOS; cols=[:Varname, :Season], rows=:Experiment,
#            sort=false))
#end
