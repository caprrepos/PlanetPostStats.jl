# Reading monthly scaleacc output files of one run into memory.

"""
    GridInfo

Horizontal and vertical grid of a run, read from the first available input files.

# Fields

  - `lon::Vector{Float32}`, `lat::Vector{Float32}`: cell-center coordinates (degrees).
  - `level::Union{Nothing,Vector{Float32}}`: layer index coordinate (`nothing` if no 3D
    variable is processed).
  - `plm`, `ple::Union{Nothing,Vector{Float32}}`: reference mid-layer and upper-edge
    pressures (hPa) of the layers, if present in the input.
  - `level_dimname::String`: name of the vertical dimension in the input (e.g. `"level"`).
  - `area::Matrix{Float64}`: grid-cell areas (`axyp`), size `(lon, lat)`.
  - `area_attrib::OrderedDict{String,Any}`: attributes of the area variable.

# Example

```julia
rd = load_run(cfg, "RUN1")
size(rd.grid.area) == (length(rd.grid.lon), length(rd.grid.lat))  # true
```
"""
struct GridInfo
    lon::Vector{Float32}
    lat::Vector{Float32}
    level::Union{Nothing,Vector{Float32}}
    plm::Union{Nothing,Vector{Float32}}
    ple::Union{Nothing,Vector{Float32}}
    level_dimname::String
    area::Matrix{Float64}
    area_attrib::OrderedDict{String,Any}
end

"""
    VariableData(data, attrib, is3d, to_pressure=true)

Monthly record of one variable of one run, held in memory.

# Fields

  - `data::Array{Float32}`: values with dimensions `(lon, lat, time)` for 2D variables or
    `(lon, lat, level, time)` for 3D variables; missing values are `NaN32`.
  - `attrib::OrderedDict{String,Any}`: output attributes (`units`, `long_name`, …).
  - `is3d::Bool`: `true` for variables with a vertical dimension.
  - `to_pressure::Bool`: for a 3D variable, whether it is interpolated to pressure levels
    when `TimeSeriesConfig.pressure_levels` is set (see [`InputVariable`](@ref)).

See also [`LazyVariable`](@ref) for variables computed only when they are processed, and
[`field`](@ref).

# Example

```julia
v = load_run(cfg, "RUN1").variables["prsurf"]
v.is3d, size(v.data)   # (false, (72, 46, 240))

w = VariableData(zeros(Float32, 4, 3, 12), OrderedDict{String,Any}("units" => "1"), false)
w.to_pressure          # true
```
"""
struct VariableData
    data::Array{Float32}
    attrib::OrderedDict{String,Any}
    is3d::Bool
    to_pressure::Bool
end

VariableData(data, attrib, is3d::Bool) = VariableData(data, attrib, is3d, true)

"""
    LazyVariable(compute, attrib, is3d, to_pressure=true)

Monthly record of one variable of one run that is computed only when it is needed:
`compute()` (a function without arguments) returns the values with the dimensions of
[`VariableData`](@ref)`.data`. The output pipeline ([`process_run`](@ref)) calls it once,
when the variable is processed, and releases the result before the next variable, so
that a run can provide many more variables than fit into memory at the same time (e.g.
optical depths for every size bin and wavelength derived from a few mass fields).

# Fields

  - `compute::Function`: `() -> AbstractArray{<:Real}`; missing values `NaN`.
  - `attrib::OrderedDict{String,Any}`: output attributes (`units`, `long_name`, …).
  - `is3d::Bool`: `true` for variables with a vertical dimension.
  - `to_pressure::Bool`: as for [`VariableData`](@ref).

# Example

```julia
mass = rand(Float32, 4, 3, 2, 12)
tau = LazyVariable(() -> 0.5f0 .* mass, OrderedDict{String,Any}("units" => "1"), true, false)
size(field(tau))       # (4, 3, 2, 12)
```
"""
struct LazyVariable
    compute::Function
    attrib::OrderedDict{String,Any}
    is3d::Bool
    to_pressure::Bool
end

LazyVariable(compute, attrib, is3d::Bool) = LazyVariable(compute, attrib, is3d, true)

"""
    AnyVariable

`Union{VariableData,LazyVariable}`: element type of the variable collections of
[`RunData`](@ref).
"""
const AnyVariable = Union{VariableData,LazyVariable}

"""
    field(v::VariableData) -> Array{Float32}
    field(v::LazyVariable) -> Array{Float32}

Values of a variable: the stored array of a [`VariableData`](@ref) (no copy), or the
result of `v.compute()` of a [`LazyVariable`](@ref), converted to `Array{Float32}`.

# Example

```julia
x = field(rd.variables["Clay"])
size(x)                # (72, 46, 40, 240)
```
"""
field(v::VariableData) = v.data
field(v::LazyVariable) = convert(Array{Float32}, v.compute())

"""
    RunData

All input data of one run, as returned by [`load_run`](@ref).

# Fields

  - `run::String`: run identifier.
  - `months::Vector{YearMonth}`: monthly record, from the first month of
    `cfg.sim_period` to the last month for which at least one input file exists.
  - `grid::GridInfo`: grid description.
  - `variables::OrderedDict{String,AnyVariable}`: variables written to the output, in
    output order ([`VariableData`](@ref) or [`LazyVariable`](@ref)).
  - `available::Matrix{Bool}`: `available[g, t]` is `true` if the input file `g` (input
    group of [`load_run`](@ref)) exists for, or contains, month `t` of `months`.
  - `auxiliary::OrderedDict{String,AnyVariable}`: variables needed for the processing
    but not written, e.g. the surface pressure for the interpolation to pressure levels
    when it is not an output variable. Variables are looked up in `variables` first.

The constructor with five arguments sets `auxiliary` empty.

# Example

```julia
rd = load_run(cfg, "RUN1")
collect(keys(rd.variables))
```
"""
struct RunData
    run::String
    months::Vector{YearMonth}
    grid::GridInfo
    variables::OrderedDict{String,AnyVariable}
    available::Matrix{Bool}
    auxiliary::OrderedDict{String,AnyVariable}
end

function RunData(run, months, grid, variables, available)
    return RunData(
        run,
        months,
        grid,
        OrderedDict{String,AnyVariable}(variables),
        available,
        OrderedDict{String,AnyVariable}(),
    )
end

"""
    lookup_variable(rd::RunData, name) -> Union{Nothing,AnyVariable}

Variable `name` of the run: from `rd.variables`, else from `rd.auxiliary`, else
`nothing`.

# Example

```julia
ps = lookup_variable(rd, "prsurf")
ps === nothing || size(field(ps))
```
"""
function lookup_variable(rd::RunData, name::AbstractString)
    haskey(rd.variables, name) && return rd.variables[name]
    haskey(rd.auxiliary, name) && return rd.auxiliary[name]
    return nothing
end

"""
    ATTRIBUTES_TO_KEEP

Input variable attributes copied to the output, besides `units` and `long_name`, if
present. Encoding attributes (`scale_factor`, `add_offset`, `_FillValue`,
`missing_value`) are deliberately *not* copied: the input is decoded when read, and the
output gets its own `_FillValue`.

# Example

```julia
julia> "standard_name" in ATTRIBUTES_TO_KEEP
true
```
"""
const ATTRIBUTES_TO_KEEP = ("standard_name", "description", "comment")

"""
    TIME_LIKE_DIMS

Dimension names treated as a (singleton) time dimension in input variables and dropped
when reading. `scaleacc` monthly files normally have no time dimension; this makes the
reader robust against files that carry a time dimension of length one.

# Example

```julia
julia> "time" in TIME_LIKE_DIMS
true
```
"""
const TIME_LIKE_DIMS = ("time", "t", "record")

"""
    spatial_shape(ds, name) -> (dims::Vector{String}, shape::Vector{Int})

Names and sizes of the spatial dimensions of variable `name` in the open dataset `ds`,
excluding singleton time-like dimensions (see [`TIME_LIKE_DIMS`](@ref)). Throws an
`ArgumentError` if the variable has fewer than 2 or more than 3 spatial dimensions.

# Example

```julia
NCDataset("JAN0001.aijlRUN1.nc") do ds
    spatial_shape(ds, "z")     # (["lon", "lat", "level"], [72, 46, 40])
end
```
"""
function spatial_shape(ds, name::AbstractString)
    var = ds[name]
    dn = collect(String, dimnames(var))
    sz = collect(Int, size(var))
    keep = [!(lowercase(d) in TIME_LIKE_DIMS && s == 1) for (d, s) in zip(dn, sz)]
    dn, sz = dn[keep], sz[keep]
    2 <= length(sz) <= 3 ||
        throw(ArgumentError("variable $name has unsupported spatial dimensions $dn"))
    return dn, sz
end

"""
    read_coordinate(ds, name) -> Union{Nothing,Vector{Float32}}

Reads the one-dimensional variable `name` from the open dataset `ds` as `Float32`
(missing values become `NaN32`); returns `nothing` if it does not exist.

# Example

```julia
NCDataset("JAN0001.aijlRUN1.nc") do ds
    read_coordinate(ds, "plm")
end
```
"""
function read_coordinate(ds, name::AbstractString)
    haskey(ds, name) || return nothing
    return Float32.(coalesce.(vec(Array(ds[name])), NaN32))
end

"""
    attribute_dict(attrib) -> OrderedDict{String,Any}

Copies a NetCDF attribute collection into a plain `OrderedDict`, detached from the
(possibly closed) file, so that it stays valid after the file is closed.

# Example

```julia
NCDataset("JAN0001.aijRUN1.nc") do ds
    attribute_dict(ds["axyp"].attrib)
end
```
"""
attribute_dict(attrib) = OrderedDict{String,Any}(String(k) => v for (k, v) in attrib)

"""
    load_run(cfg::TimeSeriesConfig, run::AbstractString) -> RunData

Reads all variables of all input groups of `run` for every month of `cfg.sim_period`
into memory, applies the variable scale factors, and evaluates the derived variables.

# Procedure

 1. Constructs the expected input path of every (group, month) and checks which files
    exist. The monthly record ends with the last month for which *any* file exists;
    months with missing files inside the record are kept as `NaN32` (gaps).
 2. Probes the first available file of each group for the variables, their
    dimensionality and attributes, and the grid (coordinates, cell area `axyp`).
 3. Reads every file once (one open per file, all variables of the group) and stores
    each variable at a fixed position, independent of which files exist.
 4. Evaluates `cfg.derived` in order.

Memory: `4 × nlon × nlat × (n2d + nlev × n3d) × nmonths` bytes for one run (e.g.
≈ 0.9 GB for 72 × 46 × 40, seven 3D variables, 240 months). Only one run is held in
memory at a time by [`run_timeseries_statistics`](@ref).

# Arguments

  - `cfg::TimeSeriesConfig`: the configuration.
  - `run::AbstractString`: run identifier (one of `cfg.runs`).

# Returns

A [`RunData`](@ref).

# Example

```julia
rd = load_run(cfg, "MarsDustRun1")
rd.months[end]                          # last month with data
size(rd.variables["Clay"].data)         # (72, 46, 40, 240)
```
"""
function load_run(cfg::TimeSeriesConfig, run::AbstractString)
    run = String(run)
    months_all = yearmonth_range(cfg.sim_period...)
    indir = joinpath(cfg.simpath, run, cfg.postdir)
    ngroups = length(cfg.groups)
    paths = [
        joinpath(indir, input_filename(cfg, ym, g.product, run)) for
        g in cfg.groups, ym in months_all
    ]
    available_all = isfile.(paths)
    any_file = vec(any(available_all; dims=1))
    last_t = findlast(any_file)
    last_t === nothing &&
        throw(ArgumentError("no input files found for run $run in $indir"))
    months = months_all[1:last_t]
    available = available_all[:, 1:last_t]
    nmissing = count(!, available)
    nmissing > 0 && @warn "Run $run: $nmissing of $(length(available)) input files " *
        "missing inside the record $(first(months)) – $(last(months)); filled with NaN."

    # --- probe the first available file of every group -------------------------------
    lon = Float32[]
    lat = Float32[]
    level = nothing
    plm = nothing
    ple = nothing
    level_dimname = "level"
    area = nothing
    area_attrib = OrderedDict{String,Any}()
    specs = Vector{Vector{Tuple{InputVariable,Vector{Int},OrderedDict{String,Any}}}}(
        undef, ngroups
    )
    for (g, group) in enumerate(cfg.groups)
        t0 = findfirst(view(available, g, :))
        t0 === nothing && throw(
            ArgumentError("run $run: no file of product $(group.product) found in $indir"),
        )
        specs[g] = NCDataset(paths[g, t0], "r") do ds
            if isempty(lon)
                lon = something(read_coordinate(ds, "lon"), Float32[])
                lat = something(read_coordinate(ds, "lat"), Float32[])
                (isempty(lon) || isempty(lat)) &&
                    throw(ArgumentError("no lon/lat coordinates in $(paths[g, t0])"))
            end
            if area === nothing && haskey(ds, cfg.area_variable)
                area = Float64.(coalesce.(Array(ds[cfg.area_variable]), 0.0))
                area_attrib = attribute_dict(ds[cfg.area_variable].attrib)
            end
            gspecs = Tuple{InputVariable,Vector{Int},OrderedDict{String,Any}}[]
            for v in group.variables
                haskey(ds, v.name) ||
                    throw(ArgumentError("variable $(v.name) not found in $(paths[g, t0])"))
                dn, shape = spatial_shape(ds, v.name)
                (shape[1], shape[2]) == (length(lon), length(lat)) || throw(
                    DimensionMismatch(
                        "variable $(v.name) has horizontal size " *
                        "$(shape[1:2]), grid is $((length(lon), length(lat)))",
                    ),
                )
                if length(shape) == 3 && level === nothing
                    level_dimname = dn[3]
                    level = something(read_coordinate(ds, dn[3]), Float32.(1:shape[3]))
                    plm = read_coordinate(ds, "plm")
                    ple = read_coordinate(ds, "ple")
                elseif length(shape) == 3 && length(level) != shape[3]
                    throw(
                        DimensionMismatch(
                            "variable $(v.name) has $(shape[3]) levels, " *
                            "expected $(length(level))",
                        ),
                    )
                end
                fileattr = ds[v.name].attrib
                att = OrderedDict{String,Any}()
                att["units"] = if !isempty(v.units)
                    v.units
                else
                    String(get(fileattr, "units", ""))
                end
                att["long_name"] = if !isempty(v.long_name)
                    v.long_name
                else
                    String(get(fileattr, "long_name", v.name))
                end
                for key in ATTRIBUTES_TO_KEEP
                    haskey(fileattr, key) && (att[key] = fileattr[key])
                end
                att["source_product"] = group.product
                push!(gspecs, (v, shape, att))
            end
            return gspecs
        end
    end
    area === nothing && throw(
        ArgumentError("area variable $(cfg.area_variable) not found in any input file")
    )
    size(area) == (length(lon), length(lat)) ||
        throw(DimensionMismatch("area variable has size $(size(area))"))
    grid = GridInfo(lon, lat, level, plm, ple, level_dimname, area, area_attrib)

    # --- allocate and read ---------------------------------------------------------------
    variables = OrderedDict{String,VariableData}()
    for g in 1:ngroups, (v, shape, att) in specs[g]
        data = fill(NaN32, shape..., last_t)
        variables[v.name] = VariableData(data, att, length(shape) == 3, v.to_pressure)
    end
    for t in 1:last_t, (g, group) in enumerate(cfg.groups)
        available[g, t] || continue
        NCDataset(paths[g, t], "r") do ds
            for (v, shape, _) in specs[g]
                haskey(ds, v.name) ||
                    throw(ArgumentError("variable $(v.name) missing in $(paths[g, t])"))
                raw = Array(ds[v.name])
                length(raw) == prod(shape) || throw(
                    DimensionMismatch(
                        "variable $(v.name) in $(paths[g, t]) has size " *
                        "$(size(raw)), expected $shape",
                    ),
                )
                dest = selectdim(variables[v.name].data, length(shape) + 1, t)
                scale = Float32(v.scale)
                dest .= coalesce.(reshape(raw, shape...), NaN32) .* scale
            end
        end
    end

    # --- derived variables ---------------------------------------------------------------
    for d in cfg.derived
        srcs = [variables[s] for s in d.sources]
        all(s -> size(s.data) == size(srcs[1].data), srcs) ||
            throw(DimensionMismatch("sources of derived variable $(d.name) differ in size"))
        data = convert(Array{Float32}, d.combine((s.data for s in srcs)...))
        size(data) == size(srcs[1].data) ||
            throw(DimensionMismatch("combine function of $(d.name) returned wrong size"))
        att = OrderedDict{String,Any}(
            "units" => isempty(d.units) ? srcs[1].attrib["units"] : d.units,
            "long_name" => d.long_name,
            "derived_from" => join(d.sources, ", "),
        )
        variables[d.name] = VariableData(data, att, srcs[1].is3d, d.to_pressure)
    end
    return RunData(run, months, grid, variables, available)
end
