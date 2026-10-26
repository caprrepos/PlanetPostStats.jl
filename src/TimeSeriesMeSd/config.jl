# Configuration types: everything that varies from simulation to simulation or from user
# to user is passed through these types instead of being hard-coded as global constants.

"""
    InputVariable(name; scale=1.0, units="", long_name="")

A variable to be read from a `scaleacc` output product (see [`InputGroup`](@ref)).

`scaleacc` stores many quantities in scaled units (e.g. tracer mixing ratios in units of
1e-8 kg kg⁻¹). After reading, values are multiplied by `scale`, so that the internal and
the output data are in physical units `units`. The output files carry no `scale_factor`
attribute for this conversion; the attribute `input_scale_factor_applied` documents it.

Any `scale_factor`/`add_offset`/`_FillValue`/`missing_value` attributes of the *input*
file are applied by NCDatasets when reading (CF decoding) before `scale` is applied.

# Arguments

  - `name::AbstractString`: variable name in the `scaleacc` file (e.g. `"prsurf"`).

# Keyword arguments

  - `scale::Real = 1.0`: multiplicative factor converting file values to `units`.
  - `units::AbstractString = ""`: units after scaling; empty → taken from the file's
    `units` attribute (only sensible if `scale == 1`).
  - `long_name::AbstractString = ""`: descriptive name; empty → the file's `long_name`.

# Example

```julia
julia> v = InputVariable("Clay"; scale=1e-8, units="kg kg-1");

julia> v.scale
1.0e-8
```
"""
struct InputVariable
    name::String
    scale::Float64
    units::String
    long_name::String
end

function InputVariable(
    name::AbstractString;
    scale::Real=1.0,
    units::AbstractString="",
    long_name::AbstractString="",
)
    isfinite(scale) && scale != 0 ||
        throw(ArgumentError("scale of $name must be finite and nonzero"))
    return InputVariable(String(name), Float64(scale), String(units), String(long_name))
end

"""
    InputGroup(product::AbstractString, variables::AbstractVector{InputVariable})

A `scaleacc` output product (e.g. `"aij"`, `"aijl"`, `"taijl"`) and the variables to be
read from its monthly files. One monthly file per product is expected for every month
of the simulation period; the file name is constructed by
`TimeSeriesConfig.filename` (default [`scaleacc_filename`](@ref)).

Variables may be two-dimensional (`lon × lat`) or three-dimensional (`lon × lat × level`);
the dimensionality is detected from the first available file.

# Example

```julia
julia> g = InputGroup("aijl", [InputVariable("z"; units="m"), InputVariable("airmass")]);

julia> [v.name for v in g.variables]
2-element Vector{String}:
 "z"
 "airmass"
```
"""
struct InputGroup
    product::String
    variables::Vector{InputVariable}
end

function InputGroup(product::AbstractString, variables::AbstractVector)
    return InputGroup(String(product), collect(InputVariable, variables))
end

"""
    DerivedVariable(name, sources, combine; units, long_name)

A variable computed from already loaded variables (input or earlier derived variables)
after reading, e.g. the accumulated dust mass mixing ratio over several size bins.

`combine` is called as `combine(a₁, a₂, …)` with the arrays of the variables named in
`sources` (in that order; `Float32`, dimensions `lon × lat [× level] × time`, missing
values `NaN32`) and must return an array of the same size. Derived variables are
evaluated in the order of `TimeSeriesConfig.derived`, so a derived variable may use
previously derived ones. All sources must have the same dimensionality.

Two helper constructors for `combine` are provided: [`sum_of`](@ref) and, for a whole
chain of accumulated sums, [`cumulative_sums`](@ref).

# Arguments

  - `name::AbstractString`: output variable name (must be unique).
  - `sources::AbstractVector{<:AbstractString}`: names of the variables combined.
  - `combine`: function of `length(sources)` arrays.

# Keyword arguments

  - `units::AbstractString = ""`: units of the result; empty → units of the first source.
  - `long_name::AbstractString = name`: descriptive name.

# Example

```julia
julia> d = DerivedVariable("ClaySilt1", ["Clay", "Silt1"], sum_of(); units="kg kg-1",
                           long_name="Dust mass mixing ratio, clay + silt bin 1");

julia> d.combine(fill(1f0, 2, 2), fill(2f0, 2, 2))
2×2 Matrix{Float32}:
 3.0  3.0
 3.0  3.0
```
"""
struct DerivedVariable{F}
    name::String
    sources::Vector{String}
    combine::F
    units::String
    long_name::String
end

function DerivedVariable(
    name::AbstractString,
    sources::AbstractVector{<:AbstractString},
    combine;
    units::AbstractString="",
    long_name::AbstractString=String(name),
)
    isempty(sources) && throw(ArgumentError("derived variable $name has no sources"))
    return DerivedVariable(
        String(name), String.(sources), combine, String(units), String(long_name)
    )
end

"""
    sum_of() -> Function

Returns a function that adds an arbitrary number of equally sized arrays element-wise
(`(a, b, …) -> a .+ b .+ …`). A `NaN` in any summand makes the sum `NaN` (missing),
i.e. a sum over size bins is only valid where all bins are valid.

# Example

```julia
julia> f = sum_of();

julia> f([1f0, 2f0], [10f0, NaN32])
2-element Vector{Float32}:
  11.0
 NaN
```
"""
sum_of() = (arrays...) -> reduce((a, b) -> a .+ b, arrays)

"""
    cumulative_sums(sources; names, long_names=names, units="") -> Vector{DerivedVariable}

Convenience constructor of a chain of accumulated sums, e.g. over size bins: for
`sources = [s₁, s₂, …, sₙ]` it returns `n - 1` derived variables

    names[1] = s₁ + s₂,  names[2] = s₁ + s₂ + s₃,  …,  names[n-1] = s₁ + … + sₙ.

Each derived variable is defined as the previous one plus the next bin, so that the
total work is linear in `n`.

# Arguments

  - `sources::AbstractVector{<:AbstractString}`: names of the size-bin variables, ordered
    from the smallest to the largest bin (`length(sources) ≥ 2`).

# Keyword arguments

  - `names`: `n - 1` names of the accumulated variables.
  - `long_names = names`: `n - 1` descriptive names.
  - `units = ""`: units (empty → units of `sources[1]`).

# Example

```julia
julia> ds = cumulative_sums(["Clay", "Silt1", "Silt2"];
                            names=["Clay_Silt1", "Clay_Silt2"], units="kg kg-1");

julia> [(d.name, d.sources) for d in ds]
2-element Vector{Tuple{String, Vector{String}}}:
 ("Clay_Silt1", ["Clay", "Silt1"])
 ("Clay_Silt2", ["Clay_Silt1", "Silt2"])
```
"""
function cumulative_sums(
    sources::AbstractVector{<:AbstractString};
    names::AbstractVector{<:AbstractString},
    long_names::AbstractVector{<:AbstractString}=names,
    units::AbstractString="",
)
    n = length(sources)
    n >= 2 || throw(ArgumentError("cumulative_sums needs at least two sources"))
    length(names) == n - 1 || throw(ArgumentError("names must have $(n - 1) elements"))
    length(long_names) == n - 1 ||
        throw(ArgumentError("long_names must have $(n - 1) elements"))
    derived = DerivedVariable[]
    previous = sources[1]
    for k in 2:n
        push!(
            derived,
            DerivedVariable(
                names[k - 1],
                [previous, sources[k]],
                sum_of();
                units=units,
                long_name=long_names[k - 1],
            ),
        )
        previous = names[k - 1]
    end
    return derived
end

"""
    dust_mass_classes(available; physvar="", addlongname="", units="") -> Vector{DerivedVariable}

Cumulative dust mass size classes of ModelE/ROCKE-3D, with the bin and class names of
`PlanetTools.DerivePhysicsVariablesModelE.accumulate_dust_masses_modele`:

| class      | sum of the bins                 |
|:-----------|:--------------------------------|
| `PM2Dust`  | `Clay`                          |
| `PM4Dust`  | `Clay + Silt1`                  |
| `PM8Dust`  | `Clay + Silt1 + Silt2`          |
| `PM16Dust` | … `+ Silt3`                     |
| `PM32Dust` | … `+ Silt4`                     |
| `PM64Dust` | … `+ Silt5`                     |

The bin variables are named `<bin><physvar>` and the classes `<class><physvar>`. As in
PlanetTools, the accumulation stops at the first bin not contained in `available`, so
only complete classes are created. The long name of a class is `"<class> <addlongname>"`
(or `"<class>"`), the units are those of `Clay<physvar>` unless `units` is given.

The bins and classes are taken from PlanetTools, so both packages stay consistent.

# Arguments

  - `available`: names of the variables that will be loaded, typically
    `[v.name for g in groups for v in g.variables]`.

# Keyword arguments

  - `physvar::AbstractString = ""`: suffix of the bin names (e.g. `"_emission"`).
  - `addlongname::AbstractString = ""`: text appended to the class name in `long_name`.
  - `units::AbstractString = ""`: units of the classes (empty → units of the clay bin).

# Example

```julia
julia> ds = dust_mass_classes(["Clay", "Silt1", "Silt2", "Silt3", "Silt4"]);

julia> [(d.name, d.sources) for d in ds][1:3]
3-element Vector{Tuple{String, Vector{String}}}:
 ("PM2Dust", ["Clay"])
 ("PM4Dust", ["PM2Dust", "Silt1"])
 ("PM8Dust", ["PM4Dust", "Silt2"])
```
"""
function dust_mass_classes(
    available;
    physvar::AbstractString="",
    addlongname::AbstractString="",
    units::AbstractString="",
)
    bins = DerivePhysicsVariablesModelE.MASS_BINS
    classes = DerivePhysicsVariablesModelE.MASS_CLASSES
    derived = DerivedVariable[]
    previous = nothing
    for (bin, class) in zip(bins, classes)
        src = bin * physvar
        src in available || break
        name = class * physvar
        long = isempty(addlongname) ? class : class * " " * addlongname
        sources = previous === nothing ? [src] : [previous, src]
        push!(
            derived, DerivedVariable(name, sources, sum_of(); units=units, long_name=long)
        )
        previous = name
    end
    return derived
end

"""
    VerticalCoordinate(ls1::Integer, psf::Real, ptop::Real)

Parameters of the hybrid sigma–pressure vertical coordinate of ModelE/ROCKE-3D, written
as scalar metadata variables `LS1`, `Psf`, and `Ptop` to all files containing 3D
variables.

Below layer `ls1` the coordinate is σ with σₑ = (pₑ − p_top)/(p_sf − p_top).

# Arguments

  - `ls1`: first layer of pure pressure coordinates above the σ layers.
  - `psf`: global mean surface pressure used for the σ definition (hPa).
  - `ptop`: pressure at the bottom of layer `ls1` (hPa).

# Example

```julia
julia> vc = VerticalCoordinate(31, 6.86, 0.2);

julia> vc.ls1
31
```
"""
struct VerticalCoordinate
    ls1::Int
    psf::Float64
    ptop::Float64
end

"""
    scaleacc_filename(ym::YearMonth, product, run; ext=".nc") -> String

Default name of a monthly `scaleacc` output file, `<MON><YYYY>.<product><run><ext>`
(e.g. `JAN0001.aijMarsDustRun1.nc`), built with
`PlanetTools.Constructors.fnames_modele` (upper-case English month abbreviation, year
zero-padded to four digits, year ≥ 0).

If your files follow a different convention, pass any function with the call signature
`(ym::YearMonth, product::String, run::String) -> String` as `TimeSeriesConfig.filename`.

# Example

```julia
julia> scaleacc_filename(YearMonth(1, 2), "taijl", "RUN1")
"FEB0001.taijlRUN1.nc"
```
"""
function scaleacc_filename(
    ym::YearMonth, product::AbstractString, run::AbstractString; ext::AbstractString=".nc"
)
    return only(
        fnames_modele(ym.year, ym.month, ym.year, ym.month, product, [run]; ext=ext)
    )
end

"""
    ALL_PRODUCTS

Tuple of all output products of [`run_timeseries_statistics`](@ref):

| symbol | content                                                     | subdirectory   |
|:-------|:------------------------------------------------------------|:---------------|
| `:mts` | monthly mean time series                                    | `MonthlyTs`    |
| `:ats` | annual mean time series (full calendar years)               | `AnnualTs`     |
| `:sts` | seasonal mean time series (one file per season)             | `SeasonalTs`   |
| `:moa` | multi-year mean and standard deviation of annual means      | `AnnualMeSd`   |
| `:mos` | multi-year mean and standard deviation of seasonal means    | `SeasonalMeSd` |
| `:mom` | multi-year mean and standard deviation of each month        | `MonthlyMeSd`  |

`:moa`, `:mos`, `:mom` each produce two files per representation, a mean (`moa`, `mos`,
`mom`) and a standard deviation file (`soa`, `sos`, `som`).

# Example

```julia
julia> :mos in ALL_PRODUCTS
true
```
"""
const ALL_PRODUCTS = (:mts, :ats, :sts, :moa, :mos, :mom)

"""
    TimeSeriesConfig(; kwargs...)

Complete, user- and simulation-specific configuration of the time-series and statistics
post-processing performed by [`run_timeseries_statistics`](@ref).

# Required keyword arguments

  - `simpath::AbstractString`: directory containing one subdirectory per run.
  - `runs::AbstractVector{<:AbstractString}`: run identifiers; input files of run `r` are
    expected in `joinpath(simpath, r, postdir)`.
  - `groups::AbstractVector{InputGroup}`: `scaleacc` products and variables to read.
  - `sim_period::Tuple{YearMonth,YearMonth}`: first and last month of the input files.
    Runs may end earlier; the record of each run ends with the last month for which at
    least one input file exists.
  - `stat_period::Tuple{YearMonth,YearMonth}`: window for the multi-year statistics.
    Only periods (years, seasons) lying entirely inside the window are used.
  - `output_dir::AbstractString`: name of the output directory created in every run
    directory (e.g. `"Tracermass3d"`).
  - `output_name::AbstractString`: tag used in output file names (e.g. `"tracermass3d"`).

# Optional keyword arguments

  - `postdir = "Post"`: subdirectory of a run directory with the `scaleacc` files.
  - `output_root = nothing`: if given, output goes to `joinpath(output_root, run,
    output_dir)` instead of `joinpath(simpath, run, output_dir)`.
  - `derived = DerivedVariable[]`: variables computed after reading (see
    [`DerivedVariable`](@ref), [`cumulative_sums`](@ref)).
  - `calendar = EARTH_NOLEAP_CALENDAR`: [`PlanetCalendar`](@ref).
  - `seasons = EARTH_SEASONS`: vector of [`SeasonSpec`](@ref).
  - `time_axis = TimeAxis()`: numerical time coordinate, see [`TimeAxis`](@ref).
  - `vertical = nothing`: [`VerticalCoordinate`](@ref) metadata for 3D output.
  - `threshold = 80.0`: validity threshold in percent (see module documentation).
  - `weight_seasons = true`: weight months by their length in seasonal means, as in
    annual means (on Mars the month lengths range from 46 to 66 sols); `false` gives
    equal weights to the months of a season.
  - `products = ALL_PRODUCTS`: subset of [`ALL_PRODUCTS`](@ref) to compute and write.
  - `fill_value = -1f30`: value written for missing data (`_FillValue`).
  - `area_variable = "axyp"`: name of the grid-cell area variable in the input files
    (the output files always call it `axyp`).
  - `global_attributes = OrderedDict{String,Any}()`: global attributes of all output
    files (the attributes `Experiment` and `File_creation_time` are added per file).
  - `filename = nothing`: function `(ym, product, run) -> String` giving the input file
    name; `nothing` → [`scaleacc_filename`](@ref).
  - `deflatelevel = 0`: NetCDF-4 compression level (0 … 9) of the output variables.
  - `overwrite = true`: overwrite existing output files; if `false`, an existing file
    raises an error.

# Example

```julia
cfg = TimeSeriesConfig(;
    simpath="/path/to/simulations",
    runs=["MarsDustRun1"],
    groups=[InputGroup("aij", [InputVariable("prsurf"; units="hPa")])],
    sim_period=(YearMonth(1, 1), YearMonth(20, 12)),
    stat_period=(YearMonth(3, 1), YearMonth(20, 12)),
    calendar=MARS_CALENDAR,
    seasons=MARS_SEASONS,
    output_dir="Prsurf",
    output_name="prsurf",
)
validate_config(cfg)
```
"""
Base.@kwdef struct TimeSeriesConfig
    simpath::String
    runs::Vector{String}
    groups::Vector{InputGroup}
    sim_period::Tuple{YearMonth,YearMonth}
    stat_period::Tuple{YearMonth,YearMonth}
    output_dir::String
    output_name::String
    postdir::String = "Post"
    output_root::Union{Nothing,String} = nothing
    derived::Vector{DerivedVariable} = DerivedVariable[]
    calendar::PlanetCalendar = EARTH_NOLEAP_CALENDAR
    seasons::Vector{SeasonSpec} = EARTH_SEASONS
    time_axis::TimeAxis = TimeAxis()
    vertical::Union{Nothing,VerticalCoordinate} = nothing
    threshold::Float64 = 80.0
    weight_seasons::Bool = true
    products::Vector{Symbol} = collect(ALL_PRODUCTS)
    fill_value::Float32 = -1.0f30
    area_variable::String = "axyp"
    global_attributes::OrderedDict{String,Any} = OrderedDict{String,Any}()
    filename::Any = nothing
    deflatelevel::Int = 0
    overwrite::Bool = true
end

"""
    input_filename(cfg::TimeSeriesConfig, ym, product, run) -> String

File name of the `scaleacc` file of `product` for month `ym` and run `run`, using
`cfg.filename` or, if that is `nothing`, [`scaleacc_filename`](@ref).

# Example

```julia
julia> input_filename(cfg, YearMonth(1, 1), "aij", "RUN1")
"JAN0001.aijRUN1.nc"
```
"""
function input_filename(cfg::TimeSeriesConfig, ym::YearMonth, product, run)
    cfg.filename === nothing && return scaleacc_filename(ym, product, run)
    return cfg.filename(ym, String(product), String(run))
end

"""
    output_directory(cfg::TimeSeriesConfig, run::AbstractString) -> String

Directory into which the output of `run` is written:
`joinpath(cfg.simpath, run, cfg.output_dir)`, or
`joinpath(cfg.output_root, run, cfg.output_dir)` if `cfg.output_root` is set.

# Example

```julia
julia> output_directory(cfg, "RUN1")
"/data/sims/RUN1/Tracermass3d"
```
"""
function output_directory(cfg::TimeSeriesConfig, run::AbstractString)
    root = cfg.output_root === nothing ? cfg.simpath : cfg.output_root
    return joinpath(root, run, cfg.output_dir)
end

"""
    variable_names(cfg::TimeSeriesConfig) -> Vector{String}

Names of all variables of the configuration, input variables first (in group order),
followed by derived variables, i.e. the order of the variables in the output files.

# Example

```julia
julia> variable_names(cfg)
3-element Vector{String}:
 "prsurf"
 "z"
 "airmass"
```
"""
function variable_names(cfg::TimeSeriesConfig)
    names = String[v.name for g in cfg.groups for v in g.variables]
    append!(names, (d.name for d in cfg.derived))
    return names
end

"""
    validate_config(cfg::TimeSeriesConfig) -> TimeSeriesConfig

Checks a configuration for consistency and throws an `ArgumentError` describing the
first problem found. Checks: non-empty runs and groups; unique variable names; sources of
derived variables defined before use; valid simulation and statistics periods, with the
statistics window overlapping the simulation period; threshold in (0, 100]; known
products; `deflatelevel` in 0:9; no reserved names (`lon`, `lat`, `level`, `time`,
`year`, `month`, `plm`, `ple`, area variable, `LS1`, `Psf`, `Ptop`) used as variable
names. Returns `cfg` so that it can be chained.

# Example

```julia
julia> validate_config(cfg);  # throws if cfg is inconsistent
```
"""
function validate_config(cfg::TimeSeriesConfig)
    isempty(cfg.runs) && throw(ArgumentError("no runs given"))
    allunique(cfg.runs) || throw(ArgumentError("run identifiers are not unique"))
    isempty(cfg.groups) && throw(ArgumentError("no input groups given"))
    allunique(g.product for g in cfg.groups) ||
        throw(ArgumentError("input groups (products) are not unique"))
    names = variable_names(cfg)
    allunique(names) || throw(ArgumentError("variable names are not unique: $names"))
    reserved = (
        "lon",
        "lat",
        "level",
        "time",
        "year",
        "month",
        "plm",
        "ple",
        cfg.area_variable,
        "axyp",
        "Dates",
        "LS1",
        "Psf",
        "Ptop",
    )
    for n in names
        n in reserved && throw(ArgumentError("variable name $n is reserved"))
    end
    known = Set(v.name for g in cfg.groups for v in g.variables)
    for d in cfg.derived
        for s in d.sources
            s in known || throw(
                ArgumentError(
                    "source $s of derived variable $(d.name) is not defined " * "before it",
                ),
            )
        end
        push!(known, d.name)
    end
    s1, s2 = cfg.sim_period
    s2 < s1 && throw(ArgumentError("sim_period ends before it starts"))
    w1, w2 = cfg.stat_period
    w2 < w1 && throw(ArgumentError("stat_period ends before it starts"))
    (w2 < s1 || s2 < w1) && throw(ArgumentError("stat_period does not overlap sim_period"))
    0 < cfg.threshold <= 100 || throw(ArgumentError("threshold must be in (0, 100]"))
    for p in cfg.products
        p in ALL_PRODUCTS || throw(ArgumentError("unknown product $p"))
    end
    0 <= cfg.deflatelevel <= 9 || throw(ArgumentError("deflatelevel must be in 0:9"))
    allunique(s.name for s in cfg.seasons) ||
        throw(ArgumentError("season names are not unique"))
    return cfg
end
