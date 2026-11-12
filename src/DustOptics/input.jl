# Reading the monthly time series of the dust tracers and building the variables of a run.

"""
    input_file(cfg::DustOpticsConfig, run) -> String

Path of the monthly time-series file with the dust tracers of `run`, written by
`TimeSeriesMeSd`: in `<input_root>/<run>/<input_dir>/MonthlyTs/` (`input_root` defaults
to `cfg.stats.simpath`), the file `<input_prefix>.<input_name><run>.<input_tag>.mts.nc`
if `cfg.input_prefix` is set, otherwise the only file whose name ends with
`.<input_name><run>.<input_tag>.mts.nc`. Throws an `ArgumentError` if no file or more
than one file is found.

# Example

```julia
julia> input_file(cfg, "MarsDustRun1")
"/path/to/simulations/MarsDustRun1/Tracermass3d/MonthlyTs/JAN0001-DEC0020.tracermass3dMarsDustRun1.ijl.mts.nc"
```
"""
function input_file(cfg::DustOpticsConfig, run::AbstractString)
    root = something(cfg.input_root, cfg.stats.simpath)
    dir = joinpath(root, run, cfg.input_dir, "MonthlyTs")
    tail = string(".", cfg.input_name, run, ".", cfg.input_tag, ".mts.nc")
    if cfg.input_prefix !== nothing
        path = joinpath(dir, cfg.input_prefix * tail)
        isfile(path) || throw(ArgumentError("input file not found: $path"))
        return path
    end
    isdir(dir) || throw(ArgumentError("run $run: input directory not found: $dir"))
    found = filter(f -> endswith(f, tail), readdir(dir))
    isempty(found) && throw(ArgumentError("run $run: no file *$tail in $dir"))
    length(found) == 1 || throw(
        ArgumentError(
            "run $run: several files *$tail in $dir: $(join(found, ", ")); " *
            "select one with input_prefix",
        ),
    )
    return joinpath(dir, only(found))
end

"""
    read_field(ds, name, tidx) -> Array{Float32}

Reads the variable `name` of the open dataset `ds` for the time indices `tidx` (last
dimension) as `Float32`; missing values (fill values) become `NaN32`.

# Example

```julia
NCDataset(path) do ds
    read_field(ds, "Clay", 1:240)          # (lon, lat, level, 240)
end
```
"""
function read_field(ds, name::AbstractString, tidx)
    haskey(ds, name) ||
        throw(ArgumentError("variable $name not found in $(NCDatasets.path(ds))"))
    var = ds[name]
    idx = ntuple(i -> i == ndims(var) ? tidx : Colon(), ndims(var))
    return Float32.(coalesce.(var[idx...], NaN32))
end

"""
    load_optics_run(cfg::DustOpticsConfig, coefficients, run) -> RunData

Reads the dust tracers and the air mass of `run` from its [`input_file`](@ref) and
returns the [`RunData`](@ref) of the optical post-processing, whose variables are
`LazyVariable`s computed when they are written (see [`DustOpticsConfig`](@ref)).

# Procedure

 1. Reads the months (`year`, `month`) of the file and keeps those within
    `cfg.stats.sim_period`, which must form a contiguous record.
 2. Reads the grid: `lon`, `lat`, the vertical coordinate of the tracers, the reference
    pressures `plm` and `ple` (if present), and the cell areas `cfg.stats.area_variable`.
 3. Reads the mass mixing ratios of the source tracers of `cfg.bins` (kg kg⁻¹) and the
    air mass `cfg.airmass` (kg m⁻²), and forms the dust mass per unit area of every
    source in every layer, `q · m_air` ([`DustBinMasses`](@ref), with the fractions of
    [`run_fractions`](@ref)). Only these source masses (and, for the density-scaled opacities, the air mass) are held in
    memory.
 4. With `cfg.stats.pressure_levels`, reads the surface pressure as auxiliary variable.
 5. For the opacities per kilometer, computes the air density of the layers from the
    layer temperature and the surface pressure ([`layer_air_density`](@ref)); the air
    mass and the air density are then also held in memory.

`coefficients` are the [`OpticalCoefficients`](@ref) of `cfg.spectral` (in order).

# Example

```julia
coefs = [interpolate_coefficients(s, cfg.bins) for s in cfg.spectral]
rd = load_optics_run(cfg, coefs, "MarsDustRun1")
collect(keys(rd.variables))[1:3]
```
"""
function load_optics_run(cfg::DustOpticsConfig, coefficients, run::AbstractString)
    run = String(run)
    path = input_file(cfg, run)
    @info "Run $run: reading $path"
    sources = bin_sources(cfg.bins)
    fractions = run_fractions(cfg, run)
    source_index = [findfirst(==(b.source), sources) for b in cfg.bins]
    pl = cfg.stats.pressure_levels
    profiles = !isempty(cfg.opacity_profiles) && !isempty(cfg.spectral)
    density_needed = needs_density(cfg)
    grid, months, masses, ps, psattr, airmass, density = NCDataset(path, "r") do ds
        for v in ("year", "month")
            haskey(ds, v) || throw(ArgumentError("variable $v not found in $path"))
        end
        allmonths = [
            YearMonth(y, m) for (y, m) in zip(Array(ds["year"]), Array(ds["month"]))
        ]
        lo, hi = cfg.stats.sim_period
        tidx = findall(m -> lo <= m <= hi, allmonths)
        isempty(tidx) &&
            throw(ArgumentError("run $run: no month of $path within the simulation period"))
        months = allmonths[tidx]
        months == yearmonth_range(first(months), last(months)) ||
            throw(ArgumentError("run $run: the months of $path are not contiguous"))
        tidx = first(tidx):last(tidx)

        haskey(ds, sources[1]) ||
            throw(ArgumentError("variable $(sources[1]) not found in $path"))
        dn = collect(String, dimnames(ds[sources[1]]))
        length(dn) == 4 || throw(
            ArgumentError(
                "$(sources[1]) in $path must have the dimensions (lon, lat, level, time)",
            ),
        )
        for v in [sources; cfg.airmass; density_needed ? cfg.temperature : String[]]
            haskey(ds, v) || throw(ArgumentError("variable $v not found in $path"))
            collect(String, dimnames(ds[v])) == dn || throw(
                DimensionMismatch(
                    "$v in $path has the dimensions $(dimnames(ds[v])), expected $dn " *
                    "(tracers, air mass, and temperature on the model levels)",
                ),
            )
        end
        lon = something(TimeSeriesMeSd.read_coordinate(ds, dn[1]), Float32[])
        lat = something(TimeSeriesMeSd.read_coordinate(ds, dn[2]), Float32[])
        nlev = size(ds[sources[1]], 3)
        level = something(TimeSeriesMeSd.read_coordinate(ds, dn[3]), Float32.(1:nlev))
        plm = TimeSeriesMeSd.read_coordinate(ds, "plm")
        ple = TimeSeriesMeSd.read_coordinate(ds, "ple")
        av = cfg.stats.area_variable
        haskey(ds, av) || throw(ArgumentError("area variable $av not found in $path"))
        area = Float64.(coalesce.(Array(ds[av]), 0.0))
        size(area) == (length(lon), length(lat)) ||
            throw(DimensionMismatch("area variable $av in $path has size $(size(area))"))
        grid = GridInfo(
            lon,
            lat,
            level,
            plm,
            ple,
            dn[3],
            area,
            TimeSeriesMeSd.attribute_dict(ds[av].attrib),
        )

        air = read_field(ds, cfg.airmass, tidx)
        masses = map(sources) do s
            q = read_field(ds, s, tidx)
            q .*= air
            return q
        end
        ps, psattr = nothing, nothing
        if pl !== nothing
            sp = pl.surface_pressure
            haskey(ds, sp) || throw(
                ArgumentError("surface pressure $sp (pressure_levels) not found in $path"),
            )
            ps = read_field(ds, sp, tidx)
            ndims(ps) == 3 || throw(ArgumentError("$sp in $path is not a 2D variable"))
            psattr = TimeSeriesMeSd.attribute_dict(ds[sp].attrib)
        end
        density = nothing
        if density_needed
            density = layer_air_density(cfg, ds, path, tidx, plm)
        end
        return grid, months, masses, ps, psattr, profiles ? air : nothing, density
    end
    m = DustBinMasses(sources, masses, source_index, fractions)
    variables = optics_variables(cfg, coefficients, m; airmass=airmass, density=density)
    auxiliary = OrderedDict{String,AnyVariable}()
    if ps !== nothing
        att = OrderedDict{String,Any}(
            k => v for (k, v) in psattr if k in ("units", "long_name")
        )
        auxiliary[pl.surface_pressure] = VariableData(ps, att, false)
    end
    available = trues(1, length(months))
    return RunData(run, months, grid, variables, available, auxiliary)
end

"""
    temperature_offset(units) -> Float64

Offset (K) that converts a temperature with the units attribute `units` to kelvin: 0 for
kelvin (`K`, `kelvin`, `degK`), 273.15 for degrees Celsius (`C`, `degC`, `deg C`, `°C`,
`degrees_C`, `degree_Celsius`, `Celsius`); case-insensitive. Throws an `ArgumentError`
for other units.

# Example

```julia
julia> temperature_offset("degC"), temperature_offset("K")
(273.15, 0.0)
```
"""
function temperature_offset(units::AbstractString)
    u = lowercase(strip(units))
    u in ("k", "kelvin", "degk", "deg k", "degrees_k") && return 0.0
    u in (
        "c",
        "degc",
        "deg c",
        "°c",
        "degrees_c",
        "degree_c",
        "degree_celsius",
        "degrees_celsius",
        "celsius",
    ) && return 273.15
    return throw(ArgumentError("unknown temperature units \"$units\" (K or °C expected)"))
end

"""
    layer_air_density(cfg::DustOpticsConfig, ds, path, tidx, plm) -> Array{Float32,4}

Air density (kg m⁻³) of the model layers of the open input dataset `ds` (file `path`)
for the time indices `tidx`, with [`air_density`](@ref) from the surface pressure
`cfg.surface_pressure` (hPa), the layer temperature `cfg.temperature` (K or °C, see
[`temperature_offset`](@ref); kelvin assumed with a warning if the variable has no
`units` attribute), the reference mid-layer pressures `plm`, `cfg.stats.vertical`, and
`cfg.gas_constant`.

# Example

```julia
NCDataset(path) do ds
    layer_air_density(cfg, ds, path, 1:240, ds["plm"][:])
end
```
"""
function layer_air_density(cfg::DustOpticsConfig, ds, path, tidx, plm)
    sp, tv = cfg.surface_pressure, cfg.temperature
    haskey(ds, sp) ||
        throw(ArgumentError("surface pressure $sp (opacity per km) not found in $path"))
    plm === nothing &&
        throw(ArgumentError("reference pressures plm (opacity per km) not found in $path"))
    ps = read_field(ds, sp, tidx)
    ndims(ps) == 3 || throw(ArgumentError("$sp in $path is not a 2D variable"))
    T = read_field(ds, tv, tidx)
    units = get(ds[tv].attrib, "units", nothing)
    offset = if units === nothing
        @warn "Temperature $tv in $path has no units attribute; kelvin assumed"
        0.0
    else
        temperature_offset(units)
    end
    iszero(offset) || (T .+= Float32(offset))
    return air_density(ps, T, plm, cfg.stats.vertical, cfg.gas_constant)
end

"""
    optics_variables(cfg::DustOpticsConfig, coefficients, m::DustBinMasses;
                     airmass=nothing, density=nothing) -> OrderedDict{String,AnyVariable}

The output variables of one run as `LazyVariable`s, in the order of
[`optics_variable_names`](@ref): each computes its field from the bin masses `m` when it
is written ([`effective_radius`](@ref), [`column_mean_effective_radius`](@ref),
[`combine_bins`](@ref) of `m` or of [`column_masses`](@ref)`(m)`,
[`opacity_per_km`](@ref), [`density_scaled_opacity`](@ref)). The column masses are
computed once, on first use. The opacity profiles need the air mass of the layers,
`airmass` (kg m⁻²), the opacities per kilometer also the air density, `density`
(kg m⁻³); both profiles are interpolated to pressure levels when
`cfg.stats.pressure_levels` is set, the optical depths of the layers are not.

# Example

```julia
vars = optics_variables(cfg, coefs, m)
size(field(vars["ext2d_550nm_PM32Dust"]))   # (lon, lat, time)
```
"""
function optics_variables(
    cfg::DustOpticsConfig, coefficients, m::DustBinMasses; airmass=nothing, density=nothing
)
    bins = cfg.bins
    colcache = Ref{Union{Nothing,DustBinMasses}}(nothing)
    function columns()
        colcache[] === nothing && (colcache[] = column_masses(m))
        return colcache[]
    end
    to2d(x) = dropdims(x; dims=3)
    vars = OrderedDict{String,AnyVariable}()
    for spec in optics_variable_specs(cfg, coefficients)
        att = OrderedDict{String,Any}("units" => spec.units, "long_name" => spec.long_name)
        k = spec.coefficients
        compute, to_p = if spec.kind == :reff3d
            (() -> effective_radius(m, bins)), cfg.effective_radius_to_pressure
        elseif spec.kind == :reffcolmass
            (() -> to2d(effective_radius(columns(), bins))), false
        elseif spec.kind == :reffcolave
            (() -> column_mean_effective_radius(m, bins)), false
        elseif spec.kind == :km
            (airmass === nothing || density === nothing) && throw(
                ArgumentError(
                    "opacity per km $(spec.name) needs the air mass and the air density",
                ),
            )
            (() -> opacity_per_km(m, k, airmass, density)), true
        elseif spec.kind == :dso
            airmass === nothing && throw(
                ArgumentError("density-scaled opacity $(spec.name) needs the air mass"),
            )
            (() -> density_scaled_opacity(m, k, airmass)), true
        elseif spec.is3d
            (() -> combine_bins(m, k)), false
        else
            (() -> to2d(combine_bins(columns(), k))), false
        end
        vars[spec.name] = LazyVariable(compute, att, spec.is3d, to_p)
    end
    return vars
end
