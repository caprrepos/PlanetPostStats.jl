# Writing the output NetCDF files.

"""
    REP_TAGS

File name tags of the spatial representations, for files with only 2D variables
(`false`) and with at least one 3D variable (`true`).

# Example

```julia
julia> REP_TAGS[true][:zm]
"zml"
```
"""
const REP_TAGS = Dict(
    false => Dict(:ij => "ij", :zm => "zm", :gm => "gm"),
    true => Dict(:ij => "ijl", :zm => "zml", :gm => "gml"),
)

"""
    REP_CELL_METHODS

Spatial part of the `cell_methods` attribute of the output variables for each spatial
representation.

# Example

```julia
julia> REP_CELL_METHODS[:gm]
"; lon, lat: mean"
```
"""
const REP_CELL_METHODS = Dict(:ij => "", :zm => "; lon: mean", :gm => "; lon, lat: mean")

"""
    OutputFile

An output NetCDF file created by [`create_output_file`](@ref). Data variables are added
with [`write_field!`](@ref), which appends to the file with
`PlanetTools.ReadWritePlanet.writenetcdf`.

# Fields

  - `path::String`: file path.
  - `rep::Symbol`: spatial representation (`:ij`, `:zm`, or `:gm`).
  - `time_cell_method::String`: temporal part of `cell_methods`
    (e.g. `"time: annual mean"`).
  - `has_time::Bool`: whether variables get a `time` dimension.
  - `level_dimname::String`: name of the vertical dimension.

# Example

```julia
of = create_output_file("x.nc", cfg, rd, :gm, true; time_cell_method="time: monthly mean",
                        time_months=rd.months)
write_field!(of, cfg, "z", gm_of_z, rd.variables["z"].attrib, true)
```
"""
struct OutputFile
    path::String
    rep::Symbol
    time_cell_method::String
    has_time::Bool
    level_dimname::String
end

"""
    software_version() -> String

Name and version of this package, written as global attribute `Software`.

# Example

```julia
julia> startswith(software_version(), "PlanetPostStats.TimeSeriesMeSd")
true
```
"""
function software_version()
    parent = parentmodule(@__MODULE__)
    v = parent === Main ? nothing : pkgversion(parent)
    return "PlanetPostStats.TimeSeriesMeSd" *
           (v === nothing ? "" : " v$v") *
           " (Julia v$VERSION)"
end

"""
    time_attributes(cfg::TimeSeriesConfig) -> OrderedDict{String,Any}

Attributes of the `time` coordinate of the time-series files: `units` from
`cfg.time_axis`, and the CF attributes of the model calendar, `calendar`
(`cfg.calendar.cf_calendar`) and `month_lengths` (CF user-defined calendar), plus
`day_length` if set. With a non-standard `calendar` (e.g. `"mars"`), CF-aware readers,
including `PlanetTools.ReadWritePlanet.readnetcdf`, return the numerical time values
instead of misinterpreting them as Earth dates.

# Example

```julia
julia> time_attributes(TimeSeriesConfig(; cfg_kwargs..., calendar=MARS_CALENDAR))["calendar"]
"mars"
```
"""
function time_attributes(cfg::TimeSeriesConfig)
    cal = cfg.calendar
    att = OrderedDict{String,Any}(
        "units" => cfg.time_axis.units,
        "long_name" => "Time (first month of averaging period)",
        "axis" => "T",
        "calendar" => cal.cf_calendar,
        "month_lengths" => collect(cal.days_in_month),
    )
    isempty(cal.day_length) || (att["day_length"] = cal.day_length)
    return att
end

"""
    create_output_file(path, cfg, rd, rep, has3d; time_cell_method, time_months=nothing,
                       climatology=false, extra_attributes=OrderedDict()) -> OutputFile

Creates the NetCDF file `path` (parent directories are created as needed) with
`PlanetTools.ReadWritePlanet.writenetcdf` and writes dimensions, coordinates, metadata,
grid-cell area, and vertical-coordinate scalars. Data variables are added afterwards with
[`write_field!`](@ref).

Contents:

  - dimensions and coordinates `lon`, `lat`, the vertical dimension (if `has3d`) and,
    for time series and climatologies, the unlimited dimension `time`;
  - for time series: `time` (months from `cfg.time_axis.reference` to the first month of
    each period, attributes from [`time_attributes`](@ref)), the integer variables
    `year(time)` and `month(time)`, and the date strings `Dates(time)`
    (`"YYYY-MM-01"`, as always written by `writenetcdf`);
  - for climatologies (`mom`, `som`): `time` = 1 … 12 and `month(time)`;
  - `plm`, `ple` (reference pressures on the vertical dimension) and `LS1`, `Psf`, `Ptop`
    (if `cfg.vertical` is set) for files with 3D variables;
  - the grid-cell area `axyp` (the name `writenetcdf` uses): cell areas (`ij`), the zonal mean of the
    cell areas (`zm`), or the total planetary area (`gm`);
  - global attributes: `Experiment`, `cfg.global_attributes`, `Planet`,
    `Missing_data_threshold_percent`, `Seasonal_means_weighted_by_month_length`,
    `Standard_deviation`, `extra_attributes`, `Software`, `File_creation_time`.

# Arguments

  - `path::AbstractString`: output file path.
  - `cfg::TimeSeriesConfig`: configuration.
  - `rd::RunData`: data of the run (grid, run identifier).
  - `rep::Symbol`: `:ij`, `:zm` (`lon` of length 1, value 0), or `:gm` (`lon` and `lat` of
    length 1, value 0).
  - `has3d::Bool`: whether the file will contain 3D variables.

# Keyword arguments

  - `time_cell_method::AbstractString`: temporal part of `cell_methods`.
  - `time_months = nothing`: first month of every time step for time-series files.
  - `climatology::Bool = false`: annual-cycle file with `time` = 1 … 12.
  - `extra_attributes`: additional global attributes.

# Example

```julia
of = create_output_file("ANN_0003_0020.prsurfRUN1.gm.moa.nc", cfg, rd, :gm, false;
                        time_cell_method="time: multi-year annual mean")
write_field!(of, cfg, "prsurf", moa, rd.variables["prsurf"].attrib, false)
```
"""
function create_output_file(
    path::AbstractString,
    cfg::TimeSeriesConfig,
    rd::RunData,
    rep::Symbol,
    has3d::Bool;
    time_cell_method::AbstractString,
    time_months::Union{Nothing,AbstractVector{YearMonth}}=nothing,
    climatology::Bool=false,
    extra_attributes=OrderedDict{String,Any}(),
)
    rep in (:ij, :zm, :gm) || throw(ArgumentError("unknown representation $rep"))
    if isfile(path)
        cfg.overwrite || throw(ArgumentError("output file exists: $path"))
        rm(path)
    end
    mkpath(dirname(path))
    grid = rd.grid
    lev = grid.level_dimname

    gattr = OrderedDict{String,Any}("Experiment" => rd.run)
    merge!(gattr, cfg.global_attributes)
    gattr["Planet"] = cfg.calendar.name
    gattr["Missing_data_threshold_percent"] = cfg.threshold
    gattr["Seasonal_means_weighted_by_month_length"] = cfg.weight_seasons ? "true" : "false"
    gattr["Standard_deviation"] = "sample (n - 1), weighted by availability"
    merge!(gattr, extra_attributes)
    gattr["Software"] = software_version()
    gattr["File_creation_time"] = string(Dates.now(UTC)) * " UTC"

    # coordinates: (name, values, attributes)
    coords = Tuple{String,Any,Any}[
        (
            "lon",
            rep == :ij ? grid.lon : Float32[0],
            OrderedDict(
                "units" => "degrees_east", "long_name" => "Longitude", "axis" => "X"
            ),
        ),
        (
            "lat",
            rep == :gm ? Float32[0] : grid.lat,
            OrderedDict(
                "units" => "degrees_north", "long_name" => "Latitude", "axis" => "Y"
            ),
        ),
    ]
    if has3d
        grid.level === nothing &&
            throw(ArgumentError("3D output requested but grid has no levels"))
        push!(coords, (lev, grid.level, OrderedDict("long_name" => "Layer", "axis" => "Z")))
    end
    has_time = time_months !== nothing || climatology
    if climatology
        tattr = OrderedDict(
            "long_name" => "Month of the annual cycle", "units" => "1", "axis" => "T"
        )
        push!(coords, ("time", Int32.(1:12), tattr))
    elseif time_months !== nothing
        tvals = Int32[time_value(cfg.time_axis, m) for m in time_months]
        push!(coords, ("time", tvals, time_attributes(cfg)))
    end
    length(coords) <= 6 || error("too many dimensions")
    ckw = Dict{Symbol,Any}()
    for (i, (_, values, attr)) in enumerate(coords)
        ckw[Symbol("datadim", i)] = values
        ckw[Symbol("dim", i, "attr")] = attr
    end

    # grid-cell area
    aattr = OrderedDict{String,Any}()
    for key in ("units", "long_name")
        haskey(grid.area_attrib, key) && (aattr[key] = grid.area_attrib[key])
    end
    areavals = if rep == :ij
        grid.area
    elseif rep == :zm
        aattr["cell_methods"] = "lon: mean"
        sum(grid.area; dims=1) ./ size(grid.area, 1)
    else
        aattr["cell_methods"] = "lon, lat: sum"
        fill(sum(grid.area), 1, 1)
    end
    datevar = if time_months === nothing || climatology
        nothing
    else
        [Date(m.year, m.month) for m in time_months]
    end
    writenetcdf(
        path;
        mode="c",
        dimnames=[c[1] for c in coords],
        ckw...,
        axypvar=Float64.(areavals),
        axypvarattr=aattr,
        datevar=datevar,
        datevarattr=OrderedDict(
            "long_name" => "Calendar date of the first day of the time period"
        ),
        globalattr=gattr,
    )

    # auxiliary time variables
    if climatology
        writenetcdf(
            path;
            mode="a",
            data=Int32.(1:12),
            names="month",
            dimnames=["time"],
            dataattr=OrderedDict("long_name" => "Month of the year"),
        )
    elseif time_months !== nothing
        writenetcdf(
            path;
            mode="a",
            data=permutedims(
                hcat(
                    Int32[m.year for m in time_months], Int32[m.month for m in time_months]
                ),
            ),
            names=["year", "month"],
            dimnames=["time"],
            dataattr=[
                OrderedDict("long_name" => "Model year of the first month of the period"),
                OrderedDict("long_name" => "First month of the period"),
            ],
        )
    end

    # vertical coordinate metadata
    if has3d
        for (name, values, long) in (
            ("plm", grid.plm, "Mid layer pressure reference level"),
            ("ple", grid.ple, "Upper edge pressure reference level"),
        )
            (values === nothing || length(values) != length(grid.level)) && continue
            writenetcdf(
                path;
                mode="a",
                data=values,
                names=name,
                dimnames=[lev],
                dataattr=OrderedDict("units" => "hPa", "long_name" => long),
            )
        end
        vc = cfg.vertical
        if vc !== nothing
            for (name, value, attr) in (
                (
                    "LS1",
                    Int32(vc.ls1),
                    OrderedDict("long_name" => "Lowest layer of pressure coordinates"),
                ),
                (
                    "Psf",
                    Float32(vc.psf),
                    OrderedDict(
                        "units" => "hPa",
                        "long_name" => "Global mean surface pressure for sigma coordinates",
                    ),
                ),
                (
                    "Ptop",
                    Float32(vc.ptop),
                    OrderedDict(
                        "units" => "hPa",
                        "long_name" => "Global mean pressure at the bottom of layer LS1",
                    ),
                ),
            )
                writenetcdf(path; mode="a", data=value, names=name, dataattr=attr)
            end
        end
    end
    return OutputFile(String(path), rep, String(time_cell_method), has_time, lev)
end

"""
    write_field!(of::OutputFile, cfg, name, data, attrib, is3d)

Appends the `Float32` variable `name` with values `data` to the output file `of` with
`PlanetTools.ReadWritePlanet.writenetcdf`.

`data` must have the dimensions `(lon, lat[, level][, time])` of the file (singleton
`lon`/`lat` for `zm`/`gm` files). Missing values (`NaN`) are written as missing, i.e. as
`cfg.fill_value`, which is set as `_FillValue`; `writenetcdf` encodes them correctly. The
attribute `cell_methods` is composed from the file's temporal method and its spatial
representation, e.g. `"time: annual mean; lon: mean"`; `attrib` (e.g. `units`,
`long_name`) is copied. The variable is compressed with `cfg.deflatelevel`.

# Example

```julia
write_field!(of, cfg, "prsurf", ats_gm, OrderedDict("units" => "hPa"), false)
```
"""
function write_field!(
    of::OutputFile,
    cfg::TimeSeriesConfig,
    name::AbstractString,
    data::AbstractArray{<:Real},
    attrib::AbstractDict,
    is3d::Bool,
)
    dims = String["lon", "lat"]
    is3d && push!(dims, of.level_dimname)
    of.has_time && push!(dims, "time")
    ndims(data) == length(dims) || throw(
        DimensionMismatch("$name: data has $(ndims(data)) dimensions, file expects $dims"),
    )
    att = OrderedDict{String,Any}(String(k) => v for (k, v) in attrib)
    att["cell_methods"] = of.time_cell_method * REP_CELL_METHODS[of.rep]
    att["_FillValue"] = cfg.fill_value
    writenetcdf(
        of.path;
        mode="a",
        data=convert(Array{Float32}, data),
        names=String(name),
        dimnames=dims,
        dataattr=att,
        deflatelevel=cfg.deflatelevel,
    )
    return nothing
end

"""
    output_filename(kind, cfg, rd, rep, has3d; prefix) -> String

Output file name `<prefix>.<output_name><run>.<tag>.<kind>.nc`, with `tag` from [`REP_TAGS`](@ref) and
`kind` one of `mts, ats, sts, moa, soa, mos, sos, mom, som`.

# Example

```julia
julia> output_filename("moa", cfg, rd, :gm, false; prefix="ANN_0003_0020")
"ANN_0003_0020.prsurfRUN1.gm.moa.nc"
```
"""
function output_filename(
    kind::AbstractString,
    cfg::TimeSeriesConfig,
    rd::RunData,
    rep::Symbol,
    has3d::Bool;
    prefix::AbstractString,
)
    return string(
        prefix, ".", cfg.output_name, rd.run, ".", REP_TAGS[has3d][rep], ".", kind, ".nc"
    )
end
