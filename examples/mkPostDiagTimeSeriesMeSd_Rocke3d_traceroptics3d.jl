# Effective radii and optical depths of the dust aerosol from ROCKE-3D Mars simulations:
# time series, multi-year means, and standard deviations.
#
# Input are the monthly time series of the 3D dust tracer mass mixing ratios and of the
# air mass of the layers written by `mkTimeSeriesMeSd_Rocke3d_tracermass3d.jl`. The clay
# tracer is split into four sub-bins with fixed mass fractions. From the dust masses of
# the size bins, the script computes the 3D effective radius of the dust, the effective
# radius of the column dust load and its mass-weighted column average, and the 3D and
# column scattering, absorption, and extinction optical depths of every bin and of the
# cumulative size classes at selected wavelengths, with the mass scattering and absorption
# coefficients of the SOCRATES optical-property files interpolated to the wavelengths and
# to the effective radii of the bins. All fields are written as time series and
# multi-year statistics of gridded, zonal-mean, and global-mean fields.
#
# Uses `PlanetPostStats.DustOptics` and `PlanetPostStats.TimeSeriesMeSd`, and
# `PlanetTools` for the size distribution of brittle fragmentation theory (add both
# packages to the active environment). All simulation- and user-specific settings are the
# constants in the section "Settings" at the top; no other code needs to be edited.
# Replace the paths, run identifiers, and metadata with your own.

using PlanetPostStats.TimeSeriesMeSd
using PlanetPostStats.DustOptics
using PlanetTools: calc_fractional_dust_masses_bft
using OrderedCollections: OrderedDict
using DataFrames

## Settings ===============================================================================

## Simulations ---------------------------------------------------------------------------
const SIMPATH = "/path/to/simulations"      # directory with one subdirectory per run
const RUNS = ["Run1", "Run2"]

## Input: output of mkTimeSeriesMeSd_Rocke3d_tracermass3d.jl ------------------------------
const GROUPDIR_IN = "Tracermass3d"          # its INPUT_DIR
const GROUPNAME_IN = "tracermass3d"         # its INPUT_NAME
const INPUT_PREFIX = nothing                # e.g. "JAN0001-DEC0020"; nothing: find file
const AIRMASS_NAME = "airmass"              # air mass of the layers (kg m-2)

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

## Optional output on pressure levels (see mkTimeSeriesMeSd_Rocke3d_tracermass3d.jl): only
## the 3D effective radius is interpolated, the optical depths of the layers are
## layer-integrated quantities. Needs the surface pressure in the input files.
const PRESSURE_LEVELS = nothing

## Statistics ----------------------------------------------------------------------------
const MISSING_THRESHOLD = 80.0   # percentage threshold required for non-missing data

## Dust size bins of the optical calculations ---------------------------------------------
const DUST_BIN_NAMES_RADIA = [
    "Clay1", "Clay2", "Clay3", "Clay4", "Silt1", "Silt2", "Silt3", "Silt4"
]
const DUST_BIN_SOURCES = [
    "Clay", "Clay", "Clay", "Clay", "Silt1", "Silt2", "Silt3", "Silt4"
]
const DUST_EFFECTIVE_RADIUS = [
    0.132e-6, 0.23e-6, 0.416e-6, 0.766e-6, 1.386e-6, 2.773e-6, 5.545e-6, 11.09e-6
]                                           # (m)
const DUST_PARTICLE_DENSITY = [2.5e3, 2.5e3, 2.5e3, 2.5e3, 2.65e3, 2.65e3, 2.65e3, 2.65e3]  # (kg m-3)

## Mass fractions of the clay tracer in the sub-bins Clay1, ..., Clay4, selected by a tag
## in the run identifier: fixed fractions, or fractions of the size distribution of
## brittle fragmentation theory (Kok, 2011) between SUBCLAY_BOUNDS.
const CLAY_SOURCE = "Clay"
const RUN_TAG_FIXED = "RAD"
const SUBCLAY_WEIGHTS_RAD = [0.009, 0.081, 0.234, 0.676]
const RUN_TAG_BFT = "AMF"
const SUBCLAY_BOUNDS = [0.1, 0.36, 0.6, 1.20, 2.00]   # particle diameters (μm)
const SUBCLAY_BFT_METHOD = :integral                   # :midpoint: earlier approximation

## Optical-property files of the radiation code and diagnostic wavelengths ---------------
const PATH_MARSDUST_SPECTRAL = "/path/to/inputfiles"
const MARSDUST_SPECTRAL_SW = "aer_sw_42_dsa_mars_sun.nc"
const MARSDUST_SPECTRAL_LW = "aer_lw_17_dsa_mars.nc"
const SPECTRAL_BANDS_SW = "sp_sw_42_dsa_mars_sun_bands.csv"
const SPECTRAL_BANDS_LW = "sp_lw_17_dsa_mars_bands.csv"
const BANDS_HEADER_ROW = 3          # line with the column names in the band files
const IDX_HUMIDITY = 1
const IDX_DUST = 1
const POST_WAVELENGTHS_SW = [0.55e-6, 0.89e-6]    # (m)
const POST_WAVELENGTHS_LW = [9.3e-6, 21.6e-6]     # (m)
const EXTRAPOLATION = :linear       # outside the tables: :linear, :flat, or :error

## Files with the coefficients interpolated to the wavelengths and bins (nothing: none)
const SPECTRAL_INTERPOL_SW = joinpath(
    PATH_MARSDUST_SPECTRAL,
    string("aer_sw_interpol_", length(POST_WAVELENGTHS_SW), "_dsa_mars_sun.nc"),
)
const SPECTRAL_INTERPOL_LW = joinpath(
    PATH_MARSDUST_SPECTRAL,
    string("aer_lw_interpol_", length(POST_WAVELENGTHS_LW), "_dsa_mars.nc"),
)

## Optical output ------------------------------------------------------------------------
const QUANTITIES = [:sct, :abs, :ext]  # scattering, absorption, extinction
const OUTPUT_LAYERS = true             # optical depths of the layers (3D)
const OUTPUT_COLUMNS = true            # column optical depths (2D)
const OUTPUT_BINS = true               # individual size bins
const OUTPUT_CLASSES = true            # cumulative size classes (SumClay1Dust, ..., PM32Dust)

## Output --------------------------------------------------------------------------------
const OUTPUT_DIR = "Traceroptics3d"     # subdirectory of each run for the output
const OUTPUT_NAME = "traceroptics3d"    # tag in the output file names

## Metadata ------------------------------------------------------------------------------
const RAW_DATA_ORIGIN = "NASA GISS ROCKE-3D Mars"
const MODEL_REPOSITORY = "<model source repository>"
const GIT_REVISION = "<git revision of the model code>"
const FILE_CREATOR = "<name>"
const ORGANIZATION = "<organization>"
const CONTACT = "<e-mail address>"

## Global means shown after the run ----------------------------------------------------
const INSPECT_VARIABLES = [
    "effRadColMass_PM32Dust",
    "effRadColAve_PM32Dust",
    "ext2d_550nm_PM32Dust",
    "abs2d_9300nm_PM32Dust",
    "ext2d_21600nm_PM32Dust",
]
const INSPECT_SEASONS = ["FMA", "ASO"]

## End of settings ========================================================================

## Size bins and mass fractions of the clay sub-bins ---------------------------------------
bins = [
    DustBin(n; source=s, radius=r, density=ρ) for (n, s, r, ρ) in zip(
        DUST_BIN_NAMES_RADIA, DUST_BIN_SOURCES, DUST_EFFECTIVE_RADIUS, DUST_PARTICLE_DENSITY
    )
]

subclay_bft = calc_fractional_dust_masses_bft(SUBCLAY_BOUNDS; method=SUBCLAY_BFT_METHOD)

function bin_fractions_of_run(run)
    weights = if occursin(RUN_TAG_FIXED, run)
        SUBCLAY_WEIGHTS_RAD
    elseif occursin(RUN_TAG_BFT, run)
        subclay_bft
    else
        error("sub-clay mass fractions not defined for run $run")
    end
    return bin_fractions(bins, CLAY_SOURCE => weights)
end

## Optical lookup tables -------------------------------------------------------------------
function spectral_set(name, ncfile, bandfile, wavelengths, coefficients_file)
    table = read_spectral_table(
        joinpath(PATH_MARSDUST_SPECTRAL, ncfile),
        joinpath(PATH_MARSDUST_SPECTRAL, bandfile);
        header_row=BANDS_HEADER_ROW,
        humidity_index=IDX_HUMIDITY,
        component_index=IDX_DUST,
    )
    return SpectralSet(name, table, wavelengths; coefficients_file=coefficients_file)
end

spectral = [
    spectral_set(
        "sw",
        MARSDUST_SPECTRAL_SW,
        SPECTRAL_BANDS_SW,
        POST_WAVELENGTHS_SW,
        SPECTRAL_INTERPOL_SW,
    ),
    spectral_set(
        "lw",
        MARSDUST_SPECTRAL_LW,
        SPECTRAL_BANDS_LW,
        POST_WAVELENGTHS_LW,
        SPECTRAL_INTERPOL_LW,
    ),
]

## Global attributes of the output files -------------------------------------------------
attributes = OrderedDict{String,Any}(
    "Raw_data_origin" => RAW_DATA_ORIGIN,
    "Model_repository" => MODEL_REPOSITORY,
    "Git_revision" => GIT_REVISION,
    "Input_created_with" => "mkTimeSeriesMeSd_Rocke3d_tracermass3d.jl",
    "Programming_language" => "Julia v" * string(VERSION),
    "Julia_script_name" => basename(@__FILE__),
    "File_creator" => FILE_CREATOR,
    "Organization" => ORGANIZATION,
    "Contact" => CONTACT,
)

## Configuration -------------------------------------------------------------------------
stats = TimeSeriesConfig(;
    simpath=SIMPATH,
    runs=RUNS,
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
)

cfg = DustOpticsConfig(;
    stats=stats,
    bins=bins,
    spectral=spectral,
    fractions=bin_fractions_of_run,
    input_dir=GROUPDIR_IN,
    input_name=GROUPNAME_IN,
    input_prefix=INPUT_PREFIX,
    airmass=AIRMASS_NAME,
    quantities=QUANTITIES,
    layers=OUTPUT_LAYERS,
    columns=OUTPUT_COLUMNS,
    bin_variables=OUTPUT_BINS,
    class_variables=OUTPUT_CLASSES,
    extrapolation=EXTRAPOLATION,
)

## Run -----------------------------------------------------------------------------------
results = run_dust_optics(cfg)

## Inspect global means --------------------------------------------------------------------
df = annual_summary(results)
filter(:Varname => in(INSPECT_VARIABLES), df)

df_seasons = seasonal_summary(results)
filter(
    [:Varname, :Season] => (v, s) -> v in INSPECT_VARIABLES && s in INSPECT_SEASONS,
    df_seasons,
)

# Optional, with SummaryTables.jl in your environment:
# using SummaryTables
# listingtable(filter(:Varname => ==("ext2d_550nm_PM32Dust"), df), :MOA; cols=[:Varname],
#     rows=:Experiment)
