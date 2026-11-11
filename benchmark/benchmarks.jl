# Benchmark suite of PlanetPostStats.jl, in the format of BenchmarkTools.jl (`SUITE`), as
# used by AirspeedVelocity.jl (`benchpkg`) and the benchmark workflow of the repository.
#
# The suite needs only BenchmarkTools and PlanetPostStats; all input data are synthetic
# and created when the script is loaded (not timed). Benchmarks of functions that do not
# exist in an older revision of the package are skipped for that revision, so that the
# same script can compare revisions.
#
# Groups:
#   TimeSeriesMeSd/kernels   statistics adapters and the pressure-level interpolation on
#                            a 72 × 46 × 40 grid, 2 Mars years of monthly fields
#   TimeSeriesMeSd/pipeline  run_timeseries_statistics: monthly scaleacc files to NetCDF
#                            time series and statistics (72 × 46 × 10 grid, 2 Mars years)
#   DustOptics/kernels       bin masses, effective radii, optical depths, opacity profiles,
#                            and the interpolation of the optical tables
#   DustOptics/pipeline      run_dust_optics on the output of the TimeSeriesMeSd pipeline
#                            (extinction of the size classes at 9.3 and 21.6 μm)
#
# Run locally (from the package root), e.g.
#   benchpkg --rev=main,dirty --bench-on=dirty
# or directly, without AirspeedVelocity:
#   julia --project=benchmark -e 'include("benchmark/benchmarks.jl"); run(SUITE; verbose=true)'

using BenchmarkTools
using PlanetPostStats
using PlanetPostStats.TimeSeriesMeSd
using PlanetPostStats.DustOptics

const NCDatasets = PlanetPostStats.TimeSeriesMeSd.NCDatasets
const OrderedDict = PlanetPostStats.TimeSeriesMeSd.OrderedDict

const SUITE = BenchmarkGroup()

# Package output (@info of the drivers) is suppressed while benchmarking.
quietly(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())

# ---------------------------------------------------------------------------------------
# Synthetic data
# ---------------------------------------------------------------------------------------

const NLON, NLAT = 72, 46
const LON = Float32.(range(-177.5, 177.5; length=NLON))
const LAT = Float32.(range(-90, 90; length=NLAT))
const AREA = [cosd(la) + 0.01 for _ in LON, la in LAT]
const MONTHS = yearmonth_range(YearMonth(1, 1), YearMonth(2, 12))
const NT = length(MONTHS)

"""
Reference mid-layer and upper-edge pressures (hPa) of `n` layers, log-spaced, and the
matching vertical coordinate: σ layers below layer `LS1 = n - n ÷ 4`, pressure layers
from there on, with `Ptop` the bottom pressure of layer `LS1`.
"""
function vertical_grid(n)
    ple = Float32.(exp.(range(log(6.5), log(0.01); length=n)))
    lower = [6.86f0; ple[1:(end - 1)]]
    plm = sqrt.(lower .* ple)
    ls1 = n - n ÷ 4
    return plm, ple, VerticalCoordinate(ls1, 6.86, ple[ls1 - 1])
end

"""Smooth positive field on (lon, lat, level, time) with some missing cells."""
function synthetic_field(nlev, scale; missing_fraction=0.002)
    x = Array{Float32}(undef, NLON, NLAT, nlev, NT)
    for t in 1:NT, l in 1:nlev, j in 1:NLAT, i in 1:NLON
        x[i, j, l, t] =
            scale * (1.5f0 + sinpi(i / 18) * cospi(j / 23) + 0.3f0 * cospi(t / 12)) /
            (1 + 0.05f0 * l)
    end
    n = round(Int, missing_fraction * length(x))
    for k in 1:n
        x[mod1(7919 * k, length(x))] = NaN32
    end
    return x
end

# ---------------------------------------------------------------------------------------
# TimeSeriesMeSd: kernels
# ---------------------------------------------------------------------------------------

let nlev = 40
    x3 = synthetic_field(nlev, 1.0f0)
    x2 = x3[:, :, 1, :]
    w = month_weights(MARS_CALENDAR, MONTHS)
    sm = planet_means(x3, AREA)
    an = annual_means(x3, sm, MONTHS, w)
    clim = climatology_slots(MONTHS, (first(MONTHS), last(MONTHS)))
    years = 1:length(annual_periods(MONTHS))
    plm, _, vc = vertical_grid(nlev)
    ps = 6.0f0 .+ x3[:, :, 1, :] ./ 10
    ps[isnan.(ps)] .= 7.0f0
    plev = default_pressure_levels("mars")

    g = SUITE["TimeSeriesMeSd"]["kernels"] = BenchmarkGroup()
    g["planet_means 2D"] = @benchmarkable planet_means($x2, $AREA)
    g["planet_means 3D"] = @benchmarkable planet_means($x3, $AREA)
    g["annual_means 3D"] = @benchmarkable annual_means($x3, $sm, $MONTHS, $w)
    g["seasonal_means 3D"] = @benchmarkable seasonal_means(
        $x3, $sm, $MONTHS, $(MARS_SEASONS[1]), $w
    )
    g["multiyear_stats 3D"] = @benchmarkable multiyear_stats($(an.ij), $(an.ij_w), $years)
    g["monthly_stats 3D"] = @benchmarkable monthly_stats($x3, nothing, $clim)
    g["interpolate_to_pressure 43 levels"] = @benchmarkable interpolate_to_pressure(
        $x3, $ps, $plm, $vc, $plev
    )
end

# ---------------------------------------------------------------------------------------
# TimeSeriesMeSd: pipeline (scaleacc files → NetCDF output)
# ---------------------------------------------------------------------------------------

const PIPE_NLEV = 10
const PIPE_ROOT = mktempdir()
const PIPE_RUN = "BenchRun"
const TRACERS = ["Clay", "Silt1", "Silt2", "Silt3", "Silt4"]

"""Writes monthly scaleacc files (aij, aijl, taijl) of one synthetic run."""
function write_scaleacc_run(root, run, nlev)
    plm, ple, _ = vertical_grid(nlev)
    temp = synthetic_field(nlev, 100.0f0; missing_fraction=0)
    air = synthetic_field(nlev, 30.0f0; missing_fraction=0)
    tracers = [
        name => synthetic_field(nlev, 200.0f0 / k) for (k, name) in enumerate(TRACERS)
    ]
    dir = joinpath(root, run, "Post")
    mkpath(dir)
    for (t, ym) in enumerate(MONTHS), product in ("aij", "aijl", "taijl")
        path = joinpath(dir, scaleacc_filename(ym, product, run))
        NCDatasets.NCDataset(path, "c") do ds
            NCDatasets.defVar(ds, "lon", LON, ("lon",))
            NCDatasets.defVar(ds, "lat", LAT, ("lat",))
            NCDatasets.defVar(ds, "axyp", AREA, ("lon", "lat"))
            if product == "aij"
                NCDatasets.defVar(
                    ds, "prsurf", 6.5f0 .+ 0.5f0 .* air[:, :, 1, t] ./ 30, ("lon", "lat")
                )
                return nothing
            end
            NCDatasets.defVar(ds, "level", Float32.(1:nlev), ("level",))
            NCDatasets.defVar(ds, "plm", plm, ("level",))
            NCDatasets.defVar(ds, "ple", ple, ("level",))
            dims = ("lon", "lat", "level")
            if product == "aijl"
                NCDatasets.defVar(ds, "TempL", 120 .+ temp[:, :, :, t], dims)
                NCDatasets.defVar(ds, "airmass", air[:, :, :, t], dims)
            else
                for (name, x) in tracers
                    NCDatasets.defVar(
                        ds,
                        name,
                        replace(x[:, :, :, t], NaN32 => -1.0f30),
                        dims;
                        attrib=["missing_value" => -1.0f30],
                    )
                end
            end
            return nothing
        end
    end
    return nothing
end

write_scaleacc_run(PIPE_ROOT, PIPE_RUN, PIPE_NLEV)

const TRACER_CFG = TimeSeriesConfig(;
    simpath=PIPE_ROOT,
    runs=[PIPE_RUN],
    groups=[
        InputGroup("aij", [InputVariable("prsurf"; units="hPa")]),
        InputGroup(
            "aijl",
            [
                InputVariable("TempL"; units="K"),
                InputVariable("airmass"; units="kg m-2", to_pressure=false),
            ],
        ),
        InputGroup(
            "taijl", [InputVariable(name; scale=1e-8, units="kg kg-1") for name in TRACERS]
        ),
    ],
    derived=dust_mass_classes(TRACERS),
    sim_period=(first(MONTHS), last(MONTHS)),
    stat_period=(first(MONTHS), last(MONTHS)),
    calendar=MARS_CALENDAR,
    seasons=MARS_SEASONS,
    vertical=vertical_grid(PIPE_NLEV)[3],
    output_dir="Tracermass3d",
    output_name="tracermass3d",
)

# The output of this run is also the input of the DustOptics pipeline.
quietly(() -> run_timeseries_statistics(TRACER_CFG))

let
    g = SUITE["TimeSeriesMeSd"]["pipeline"] = BenchmarkGroup()
    g["run_timeseries_statistics"] = @benchmarkable(
        quietly(() -> run_timeseries_statistics($TRACER_CFG)),
        evals=1,
        samples=3,
        seconds=30
    )
    cfgp = TimeSeriesConfig(;
        (
            f => getfield(TRACER_CFG, f) for
            f in fieldnames(TimeSeriesConfig) if f ∉ (:pressure_levels, :output_dir)
        )...,
        pressure_levels=PressureLevels(),
        output_dir="Tracermass3dPlev",
    )
    g["run_timeseries_statistics pressure levels"] = @benchmarkable(
        quietly(() -> run_timeseries_statistics($cfgp)), evals=1, samples=3, seconds=30
    )
end

# ---------------------------------------------------------------------------------------
# DustOptics
# ---------------------------------------------------------------------------------------

const BINS = [
    DustBin("Clay1"; source="Clay", radius=0.132e-6, density=2.5e3),
    DustBin("Clay2"; source="Clay", radius=0.23e-6, density=2.5e3),
    DustBin("Clay3"; source="Clay", radius=0.416e-6, density=2.5e3),
    DustBin("Clay4"; source="Clay", radius=0.766e-6, density=2.5e3),
    DustBin("Silt1"; radius=1.386e-6, density=2.65e3),
    DustBin("Silt2"; radius=2.773e-6, density=2.65e3),
    DustBin("Silt3"; radius=5.545e-6, density=2.65e3),
    DustBin("Silt4"; radius=11.09e-6, density=2.65e3),
]
const FRACTIONS = bin_fractions(BINS, "Clay" => [0.009, 0.081, 0.234, 0.676])

"""Optical table of 42 bands and 30 radii following power laws."""
function synthetic_table()
    edges = exp.(range(log(0.2e-6), log(50e-6); length=43))
    bounds = [(edges[b], edges[b + 1]) for b in 1:42]
    centers = sqrt.(edges[1:(end - 1)] .* edges[2:end])
    radii = exp.(range(log(0.05e-6), log(30e-6); length=30))
    tab(f) = [f(λ, r) for λ in centers, r in radii]
    return SpectralTable(
        bounds,
        radii,
        tab((λ, r) -> 2e3 * (λ / 1e-6)^-0.5 * (r / 1e-6)^-1),
        tab((λ, r) -> 3e2 * (λ / 1e-6)^0.7 * (r / 1e-6)^-0.8),
        tab((λ, r) -> 0.7 * (λ / 1e-6)^-0.1),
    )
end
const SPECTRAL = [
    SpectralSet("sw", synthetic_table(), [0.55e-6]),
    SpectralSet("lw", synthetic_table(), [9.3e-6, 21.6e-6]),
]

let nlev = 40
    sources = ["Clay", "Silt1", "Silt2", "Silt3", "Silt4"]
    air = synthetic_field(nlev, 30.0f0; missing_fraction=0)
    masses = [synthetic_field(nlev, 1.0f-6 / k) .* air for k in 1:5]
    m = DustBinMasses(sources, masses, [1, 1, 1, 1, 2, 3, 4, 5], FRACTIONS)
    k = collect(range(100.0, 800.0; length=8))
    plm, _, vc = vertical_grid(nlev)
    ps = fill(6.86f0, NLON, NLAT, NT)
    temp = 120 .+ synthetic_field(nlev, 100.0f0; missing_fraction=0)

    g = SUITE["DustOptics"]["kernels"] = BenchmarkGroup()
    g["combine_bins"] = @benchmarkable combine_bins($m, $k)
    g["column_masses"] = @benchmarkable column_masses($m)
    g["effective_radius"] = @benchmarkable effective_radius($m, $BINS)
    g["column_mean_effective_radius"] = @benchmarkable column_mean_effective_radius(
        $m, $BINS
    )
    g["interpolate_coefficients"] = @benchmarkable interpolate_coefficients(
        $(SPECTRAL[2]), $BINS
    )
    if isdefined(DustOptics, :density_scaled_opacity)
        g["density_scaled_opacity"] = @benchmarkable DustOptics.density_scaled_opacity(
            $m, $k, $air
        )
    end
    if isdefined(DustOptics, :air_density) && isdefined(DustOptics, :opacity_per_km)
        rho = DustOptics.air_density(ps, temp, plm, vc, 188.92)
        g["air_density"] = @benchmarkable DustOptics.air_density(
            $ps, $temp, $plm, $vc, 188.92
        )
        g["opacity_per_km"] = @benchmarkable DustOptics.opacity_per_km($m, $k, $air, $rho)
    end
end

let
    stats = TimeSeriesConfig(;
        simpath=PIPE_ROOT,
        runs=[PIPE_RUN],
        sim_period=(first(MONTHS), last(MONTHS)),
        stat_period=(first(MONTHS), last(MONTHS)),
        calendar=MARS_CALENDAR,
        seasons=MARS_SEASONS,
        vertical=TRACER_CFG.vertical,
        output_dir="Traceroptics3d",
        output_name="traceroptics3d",
    )
    cfg = DustOpticsConfig(;
        stats=stats,
        bins=BINS,
        spectral=[SPECTRAL[2]],
        fractions=FRACTIONS,
        quantities=[:ext],
        bin_variables=false,
        input_dir="Tracermass3d",
        input_name="tracermass3d",
    )
    g = SUITE["DustOptics"]["pipeline"] = BenchmarkGroup()
    g["run_dust_optics"] = @benchmarkable(
        quietly(() -> run_dust_optics($cfg)), evals=1, samples=3, seconds=30
    )
end
