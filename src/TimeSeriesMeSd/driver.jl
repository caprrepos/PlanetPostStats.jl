# Top-level driver: processes runs one after the other.

"""
    RunResult

Result of processing one run with [`run_timeseries_statistics`](@ref).

# Fields

  - `run::String`: run identifier.
  - `months::Vector{YearMonth}`: monthly record that was processed.
  - `files::Vector{String}`: paths of all files written.
  - `annual::DataFrame`: multi-year statistics of the annual *global* means of all 2D
    variables, columns `Varname`, `N_years`, `MOA`, `SOA`, `Experiment`, for interactive
    inspection.
  - `seasonal::DataFrame`: same for seasonal global means, columns `Varname`, `Season`,
    `N_seasons`, `MOS`, `SOS`, `Experiment`.

# Example

```julia
res = run_timeseries_statistics(cfg)
first(res).annual
```
"""
struct RunResult
    run::String
    months::Vector{YearMonth}
    files::Vector{String}
    annual::DataFrame
    seasonal::DataFrame
end

"""
    annual_summary(results::AbstractVector{RunResult}) -> DataFrame

Concatenation of the `annual` tables of all runs (multi-year mean `MOA` and standard
deviation `SOA` of the annual global means of the 2D variables).

# Example

```julia
using DataFrames
df = annual_summary(run_timeseries_statistics(cfg))
filter(:Varname => ==("prsurf"), df)
```
"""
function annual_summary(results::AbstractVector{RunResult})
    return reduce(vcat, (r.annual for r in results))
end

"""
    seasonal_summary(results::AbstractVector{RunResult}) -> DataFrame

Concatenation of the `seasonal` tables of all runs (multi-year mean `MOS` and standard
deviation `SOS` of the seasonal global means of the 2D variables).

# Example

```julia
using DataFrames
df = seasonal_summary(run_timeseries_statistics(cfg))
filter([:Varname, :Season] => (v, s) -> v == "prsurf" && s == "FMA", df)
```
"""
function seasonal_summary(results::AbstractVector{RunResult})
    return reduce(vcat, (r.seasonal for r in results))
end

"""
    period_label(months, periods) -> String

`"YYYY-YYYY"` label with the years of the first months of the first and last period.

# Example

```julia
julia> months = yearmonth_range(YearMonth(1, 1), YearMonth(3, 12));

julia> period_label(months, annual_periods(months))
"0001-0003"
```
"""
function period_label(months, periods::AbstractVector{Period}, sep::AbstractString="-")
    return string(yyyy(first(periods).first.year), sep, yyyy(last(periods).first.year))
end

"""
    run_timeseries_statistics(cfg::TimeSeriesConfig; runs=cfg.runs) -> Vector{RunResult}

Main entry point. Validates `cfg` and processes every run in `runs` (default: all runs of
the configuration) sequentially:

 1. [`load_run`](@ref) reads the monthly `scaleacc` files and evaluates derived
    variables.
 2. With `cfg.pressure_levels`, the 3D variables with `to_pressure = true` are
    interpolated to pressure levels ([`interpolate_to_pressure`](@ref)), one variable
    at a time, and assigned to the output files ([`output_entries`](@ref)).
 3. For every variable, the monthly gridded field (`ij`), its zonal mean (`zm`), and its
    global mean (`gm`) are formed ([`planet_means`](@ref)); for each
    of the three, the requested products of [`ALL_PRODUCTS`](@ref) are computed
    ([`annual_means`](@ref), [`seasonal_means`](@ref), [`multiyear_stats`](@ref),
    [`monthly_stats`](@ref), all built on `PlanetTools.PlanetStatistics`)
    and written immediately, so that only one variable's results are held in memory.
 4. The run's data are released before the next run is loaded.

Output directory of run `r`: [`output_directory`](@ref)`(cfg, r)`, with subdirectories
`MonthlyTs`, `AnnualTs`, `SeasonalTs`, `AnnualMeSd`, `SeasonalMeSd`, `MonthlyMeSd`.

Statistics are computed only from periods lying entirely in `cfg.stat_period`; products
for which no such period exists in a run (e.g. a run shorter than the window) are
skipped with a warning.

# Arguments

  - `cfg::TimeSeriesConfig`: configuration.

# Keyword arguments

  - `runs = cfg.runs`: subset of runs to process.

# Returns

`Vector{RunResult}`, one per run, with written file paths and summary tables
(see [`annual_summary`](@ref), [`seasonal_summary`](@ref)).

# Example

```julia
using PlanetPostStats.TimeSeriesMeSd
results = run_timeseries_statistics(cfg)
annual_summary(results)
```
"""
function run_timeseries_statistics(cfg::TimeSeriesConfig; runs=cfg.runs)
    validate_config(cfg)
    return [process_run(cfg, String(r)) for r in runs]
end

"""
    OutputEntry

One variable as written to one set of output files.

# Fields

  - `fileset::Symbol`: `:main` (regular files) or `:pressure` (additional files on
    pressure levels, `PressureLevels(…; mode=:add)`).
  - `name::String`: variable name.
  - `attrib::OrderedDict{String,Any}`: output attributes.
  - `vertical::Symbol`: `:none` (2D), `:model` (model levels), or `:pressure`.
  - `data::Function`: function without arguments returning the monthly field; for
    pressure levels it performs the interpolation, so that interpolated fields exist only
    while the variable is processed.

# Example

```julia
entries = output_entries(cfg, load_run(cfg, "MarsDustRun1"))
[(e.fileset, e.name, e.vertical) for e in entries]
```
"""
struct OutputEntry
    fileset::Symbol
    name::String
    attrib::OrderedDict{String,Any}
    vertical::Symbol
    data::Function
end

"""
    output_entries(cfg::TimeSeriesConfig, rd::RunData) -> Vector{OutputEntry}

Assigns the variables of a run to the output file sets and vertical coordinates:

  - without `cfg.pressure_levels`: all variables in `:main`, 3D variables on model levels;
  - `mode = :add`: all variables in `:main` as above, plus every 3D variable with
    `to_pressure = true` on pressure levels in `:pressure`;
  - `mode = :replace`: all variables in `:main`, 3D variables with `to_pressure = true` on
    pressure levels, the other 3D variables on model levels.

The pressure-level fields are computed with [`interpolate_to_pressure`](@ref) from the
monthly surface pressure `cfg.pressure_levels.surface_pressure` of the run.

# Example

```julia
entries = output_entries(cfg, rd)
count(e -> e.vertical == :pressure, entries)
```
"""
function output_entries(cfg::TimeSeriesConfig, rd::RunData)
    pl = cfg.pressure_levels
    plev = pl === nothing ? nothing : output_pressure_levels(cfg)
    flags = pressure_flags(cfg)
    main = OutputEntry[]
    extra = OutputEntry[]
    ps = nothing
    if pl !== nothing
        haskey(rd.variables, pl.surface_pressure) || throw(
            ArgumentError("surface-pressure variable $(pl.surface_pressure) not loaded")
        )
        psvar = rd.variables[pl.surface_pressure]
        psvar.is3d && throw(ArgumentError("$(pl.surface_pressure) is not a 2D variable"))
        ps = psvar.data
        rd.grid.plm === nothing &&
            throw(ArgumentError("no reference pressures plm in the input files"))
    end
    for (name, v) in rd.variables
        data = v.data
        if pl !== nothing && v.is3d && get(flags, name, true)
            interp =
                () -> interpolate_to_pressure(data, ps, rd.grid.plm, cfg.vertical, plev)
            if pl.mode == :add
                push!(main, OutputEntry(:main, name, v.attrib, :model, () -> data))
                push!(extra, OutputEntry(:pressure, name, v.attrib, :pressure, interp))
            else
                push!(main, OutputEntry(:main, name, v.attrib, :pressure, interp))
            end
        else
            vertical = v.is3d ? :model : :none
            push!(main, OutputEntry(:main, name, v.attrib, vertical, () -> data))
        end
    end
    return vcat(main, extra)
end

"""
    process_run(cfg::TimeSeriesConfig, run::String) -> RunResult

Processes a single run (see [`run_timeseries_statistics`](@ref) for the steps).

# Example

```julia
res = process_run(cfg, "MarsDustRun1")
res.files
```
"""
function process_run(cfg::TimeSeriesConfig, run::String)
    @info "Reading run $run"
    rd = load_run(cfg, run)
    months = rd.months
    nbytes = sum(sizeof(v.data) for v in values(rd.variables))
    @info "Run $run: $(length(months)) months ($(first(months)) – $(last(months))), " *
        "$(length(rd.variables)) variables, $(round(nbytes / 2^20; digits=1)) MiB"

    thr = cfg.threshold
    w = month_weights(cfg.calendar, months)
    products = Set(cfg.products)
    entries = output_entries(cfg, rd)
    window = cfg.stat_period
    cal = cfg.calendar
    outdir = output_directory(cfg, run)

    years = annual_periods(months)
    stat_years = findall(p -> period_in_window(p, months, window), years)
    seasons = [seasonal_periods(months, s) for s in cfg.seasons]
    stat_seasons = [findall(p -> period_in_window(p, months, window), ps) for ps in seasons]
    clim = climatology_slots(months, window)

    need_ats = :ats in products || :moa in products
    need_sts = :sts in products || :mos in products
    do_moa = :moa in products && !isempty(stat_years)
    do_mos = [:mos in products && !isempty(stat_seasons[k]) for k in eachindex(seasons)]
    do_mom = :mom in products && any(!isempty, clim)
    need_ats &&
        isempty(years) &&
        @warn "Run $run: no complete year, annual products skipped"
    :moa in products &&
        !isempty(years) &&
        isempty(stat_years) &&
        @warn "Run $run: no complete year inside stat_period, moa/soa skipped"
    :mom in products &&
        !do_mom &&
        @warn "Run $run: no complete year inside stat_period, mom/som skipped"
    for k in eachindex(seasons)
        need_sts &&
            isempty(seasons[k]) &&
            @warn "Run $run: no complete season $(cfg.seasons[k].name)"
        :mos in products &&
            !isempty(seasons[k]) &&
            !do_mos[k] &&
            @warn "Run $run: no season $(cfg.seasons[k].name) inside stat_period"
    end

    files = OrderedDict{Any,OutputFile}()
    sub(d) = joinpath(outdir, d)
    filesets = unique(e.fileset for e in entries)
    has_model = Dict(
        fs => any(e -> e.fileset == fs && e.vertical == :model, entries) for fs in filesets
    )
    has_plev = Dict(
        fs => any(e -> e.fileset == fs && e.vertical == :pressure, entries) for
        fs in filesets
    )
    vtag(fs) =
        if has_model[fs]
            "l"
        elseif has_plev[fs]
            "p"
        else
            ""
        end
    function open!(key, dir, kind, rep, prefix, tcm; kw...)
        fs = first(key)
        name = output_filename(kind, cfg, rd, rep, vtag(fs); prefix=prefix)
        files[key] = create_output_file(
            joinpath(sub(dir), name),
            cfg,
            rd,
            rep,
            has_model[fs];
            has_plev=has_plev[fs],
            time_cell_method=tcm,
            kw...,
        )
        return nothing
    end
    function statattr(n, what, label)
        return OrderedDict{String,Any}(
            "Statistics_window" =>
                "$(yyyy(window[1].year))-$(lpad(window[1].month, 2, '0'))" *
                " to $(yyyy(window[2].year))-$(lpad(window[2].month, 2, '0'))",
            "Statistics_periods_used" => label,
            "Number_of_$what" => n,
        )
    end

    annual = DataFrame(;
        Varname=String[], N_years=Int[], MOA=Float32[], SOA=Float32[], Experiment=String[]
    )
    seasonal = DataFrame(;
        Varname=String[],
        Season=String[],
        N_seasons=Int[],
        MOS=Float32[],
        SOS=Float32[],
        Experiment=String[],
    )

    begin
        # ---- create all output files of the run ------------------------------------------
        mprefix = string(
            monthname(cal, first(months)),
            yyyy(first(months).year),
            "-",
            monthname(cal, last(months)),
            yyyy(last(months).year),
        )
        for fs in filesets, rep in (:ij, :zm, :gm)
            if :mts in products
                open!(
                    (fs, :mts, rep),
                    "MonthlyTs",
                    "mts",
                    rep,
                    mprefix,
                    "time: monthly mean";
                    time_months=months,
                )
            end
            if :ats in products && !isempty(years)
                open!(
                    (fs, :ats, rep),
                    "AnnualTs",
                    "ats",
                    rep,
                    "ANN_" * period_label(months, years),
                    "time: annual mean";
                    time_months=[p.first for p in years],
                )
            end
            if do_moa
                label = period_label(months, years[stat_years], "_")
                att = statattr(length(stat_years), "years", label)
                open!(
                    (fs, :moa, rep),
                    "AnnualMeSd",
                    "moa",
                    rep,
                    "ANN_" * label,
                    "time: multi-year annual mean";
                    extra_attributes=att,
                )
                open!(
                    (fs, :soa, rep),
                    "AnnualMeSd",
                    "soa",
                    rep,
                    "ANN_" * label,
                    "time: multi-year annual standard deviation";
                    extra_attributes=att,
                )
            end
            for (k, s) in enumerate(cfg.seasons)
                ps = seasons[k]
                if :sts in products && !isempty(ps)
                    open!(
                        (fs, :sts, rep, k),
                        "SeasonalTs",
                        "sts",
                        rep,
                        s.name * "_" * period_label(months, ps),
                        "time: seasonal mean";
                        time_months=[p.first for p in ps],
                    )
                end
                if do_mos[k]
                    label = period_label(months, ps[stat_seasons[k]], "_")
                    att = statattr(length(stat_seasons[k]), "seasons", label)
                    open!(
                        (fs, :mos, rep, k),
                        "SeasonalMeSd",
                        "mos",
                        rep,
                        s.name * "_" * label,
                        "time: multi-year seasonal mean";
                        extra_attributes=att,
                    )
                    open!(
                        (fs, :sos, rep, k),
                        "SeasonalMeSd",
                        "sos",
                        rep,
                        s.name * "_" * label,
                        "time: multi-year seasonal standard deviation";
                        extra_attributes=att,
                    )
                end
            end
            if do_mom
                used = filter(p -> period_in_window(p, months, window), years)
                label = period_label(months, used, "_")
                prefix = cal.month_names[1] * "-" * cal.month_names[12] * "_" * label
                att = statattr(length(used), "years", label)
                open!(
                    (fs, :mom, rep),
                    "MonthlyMeSd",
                    "mom",
                    rep,
                    prefix,
                    "time: multi-year monthly mean";
                    climatology=true,
                    extra_attributes=att,
                )
                open!(
                    (fs, :som, rep),
                    "MonthlyMeSd",
                    "som",
                    rep,
                    prefix,
                    "time: multi-year monthly standard deviation";
                    climatology=true,
                    extra_attributes=att,
                )
            end
        end

        # ---- variable by variable ----------------------------------------------------------
        stat_year_range = contiguous_range(stat_years)
        stat_season_ranges = [contiguous_range(ix) for ix in stat_seasons]
        seasw = cfg.weight_seasons ? w : nothing
        for e in entries
            name = e.name
            fs = e.fileset
            @info "Run $run: processing $name" * (e.vertical == :pressure ? " (plev)" : "")
            att = e.attrib
            is3d = e.vertical != :none
            x = e.data()
            sm = planet_means(x, rd.grid.area; threshold=thr)
            # (field, availability) of the monthly ij, zm, gm representations; the gridded
            # monthly field enters with availability `nothing` (equal weights)
            monthly = (
                :ij => (x, nothing), :zm => (sm.zm, sm.zm_w), :gm => (sm.gm, sm.gm_w)
            )
            write!(key, data) =
                write_field!(files[(fs, key...)], cfg, name, data, att, e.vertical)

            if :mts in products
                for (rep, (xr, _)) in monthly
                    write!((:mts, rep), xr)
                end
            end

            if need_ats && !isempty(years)
                an = annual_means(x, sm, months, w; threshold=thr)
                annual_reps = (
                    :ij => (an.ij, an.ij_w),
                    :zm => (an.zm, an.zm_w),
                    :gm => (an.gm, an.gm_w),
                )
                for (rep, (xr, xw)) in annual_reps
                    :ats in products && write!((:ats, rep), xr)
                    do_moa || continue
                    me, sd, _ = multiyear_stats(xr, xw, stat_year_range; threshold=thr)
                    write!((:moa, rep), me)
                    write!((:soa, rep), sd)
                    if rep == :gm && !is3d
                        push!(annual, (name, length(stat_years), me[1], sd[1], run))
                    end
                end
                an = nothing
            end

            for (k, s) in enumerate(cfg.seasons)
                (need_sts && !isempty(seasons[k])) || continue
                se = seasonal_means(x, sm, months, s, seasw; threshold=thr)
                season_reps = (
                    :ij => (se.ij, se.ij_w),
                    :zm => (se.zm, se.zm_w),
                    :gm => (se.gm, se.gm_w),
                )
                for (rep, (xr, xw)) in season_reps
                    :sts in products && write!((:sts, rep, k), xr)
                    do_mos[k] || continue
                    me, sd, _ = multiyear_stats(
                        xr, xw, stat_season_ranges[k]; threshold=thr
                    )
                    write!((:mos, rep, k), me)
                    write!((:sos, rep, k), sd)
                    if rep == :gm && !is3d
                        push!(
                            seasonal,
                            (name, s.name, length(stat_seasons[k]), me[1], sd[1], run),
                        )
                    end
                end
            end

            if do_mom
                for (rep, (xr, xw)) in monthly
                    me, sd, _ = monthly_stats(xr, xw, clim; threshold=thr)
                    write!((:mom, rep), me)
                    write!((:som, rep), sd)
                end
            end
        end
    end
    paths = [of.path for of in values(files)]
    @info "Run $run: wrote $(length(paths)) files to $outdir"
    return RunResult(run, months, paths, annual, seasonal)
end
