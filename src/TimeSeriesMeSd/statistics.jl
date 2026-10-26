# Thin adapters to PlanetTools.PlanetStatistics.
#
# PlanetTools expects arrays laid out as (Sim, Var, Lon, Lat, [Lev], Time). This package
# processes one variable of one run at a time with arrays (Lon, Lat, [Lev], Time), so the
# adapters add two leading singleton dimensions with `reshape` (no copy), call the
# PlanetTools function, and drop the singleton dimensions again. All statistical
# definitions (weighting, availability, threshold, sample standard deviation) are therefore
# exactly those of PlanetTools.

"""
    lift(a) -> AbstractArray
    lift(nothing) -> nothing

Adds two leading singleton dimensions `(Sim = 1, Var = 1)` to `a` without copying, for
calls of PlanetTools functions. `nothing` is passed through.

# Example

```julia
julia> size(lift(zeros(4, 3, 10)))
(1, 1, 4, 3, 10)
```
"""
lift(a::AbstractArray) = reshape(a, 1, 1, size(a)...)
lift(::Nothing) = nothing

"""
    unlift(a, n=2) -> AbstractArray

Removes the `n` leading singleton dimensions of `a` (inverse of [`lift`](@ref); `n = 3`
also removes the `Season` dimension of seasonal arrays), without copying.

# Example

```julia
julia> size(unlift(zeros(1, 1, 4, 3)))
(4, 3)
```
"""
function unlift(a::AbstractArray, n::Integer=2)
    all(==(1), size(a)[1:n]) || throw(ArgumentError("leading dimensions are not singleton"))
    return reshape(a, size(a)[(n + 1):end])
end

"""
    nanarray(dims...) -> Array{Float32}

`Float32` array of size `dims` filled with `NaN32`. PlanetTools leaves output entries not
covered by a calculation unchanged, so outputs are initialized with `NaN`.

# Example

```julia
julia> all(isnan, nanarray(2, 3))
true
```
"""
nanarray(dims...) = fill(NaN32, dims...)

"""
    SpatialMeans

Zonal and global means of a monthly field with their availability weights, as returned
by [`planet_means`](@ref).

# Fields

  - `zm`, `zm_w`: zonal mean and fraction of valid longitudes, `(1, Lat, [Lev], Time)`.
  - `gm`, `gm_w`: area-weighted global mean and valid area fraction, `(1, 1, [Lev], Time)`.

Availabilities are `0` where the threshold is not met (the mean is then `NaN`).

# Example

```julia
sm = planet_means(x, area; threshold=80)
size(sm.gm)
```
"""
struct SpatialMeans
    zm::Array{Float32}
    zm_w::Array{Float32}
    gm::Array{Float32}
    gm_w::Array{Float32}
end

"""
    planet_means(x, area; threshold=80.0) -> SpatialMeans

Zonal means and area-weighted global means of the monthly field `x`, computed with
`PlanetTools.PlanetStatistics.calc_planet_means!`.

Definitions (PlanetTools): the zonal mean averages the valid longitudes with equal
weights, its availability is the fraction of valid longitudes; the global mean is
weighted by the cell areas `area`, its availability is the valid area fraction. A mean is
`NaN` (availability 0) if the availability is below `threshold` percent.

The PlanetTools kernel loops over the data without temporary arrays and is
multi-threaded over levels and time steps.

# Arguments

  - `x::AbstractArray{<:Real}`: `(Lon, Lat, Time)` or `(Lon, Lat, Lev, Time)`, missing
    values `NaN`.
  - `area::AbstractMatrix{<:Real}`: cell areas `(Lon, Lat)` (`axyp`).

# Keyword arguments

  - `threshold::Real = 80.0`: minimum valid percentage (longitudes or area).

# Example

```julia
julia> x = reshape(Float32[1, 3, 2, NaN], 2, 2, 1);

julia> sm = planet_means(x, [1.0 3.0; 1.0 3.0]; threshold=0);

julia> vec(sm.zm), sm.gm[1], sm.gm_w[1]
(Float32[2.0, 2.0], 2.0f0, 0.625f0)
```
"""
function planet_means(
    x::AbstractArray{<:Real}, area::AbstractMatrix{<:Real}; threshold::Real=80.0
)
    nd = ndims(x)
    3 <= nd <= 4 || throw(ArgumentError("x must be (Lon, Lat, [Lev], Time)"))
    sz = size(x)
    nt = sz[end]
    zm = nanarray(1, sz[2:end]...)
    zm_w = nanarray(1, sz[2:end]...)
    gm = nanarray(1, 1, sz[3:end]...)
    gm_w = nanarray(1, 1, sz[3:end]...)
    calc_planet_means!(
        lift(x), lift(zm), lift(gm), lift(zm_w), lift(gm_w), area, [nt]; threshold=threshold
    )
    return SpatialMeans(zm, zm_w, gm, gm_w)
end

"""
    PeriodMeans

Temporal means over periods (years or occurrences of one season) of the gridded field and
of its zonal and global means, with availability weights (fraction of valid months of
the period; `0` where the threshold is not met). Returned by [`annual_means`](@ref) and
[`seasonal_means`](@ref).

# Fields

  - `ij`, `ij_w`: `(Lon, Lat, [Lev], Period)`.
  - `zm`, `zm_w`: `(1, Lat, [Lev], Period)`.
  - `gm`, `gm_w`: `(1, 1, [Lev], Period)`.

# Example

```julia
an = annual_means(x, sm, months, w; threshold=80)
size(an.gm, ndims(an.gm))     # number of years
```
"""
struct PeriodMeans
    ij::Array{Float32}
    ij_w::Array{Float32}
    zm::Array{Float32}
    zm_w::Array{Float32}
    gm::Array{Float32}
    gm_w::Array{Float32}
end

"""
    yearmonth_dates(months::AbstractVector{YearMonth}) -> Vector{Dates.Date}

Converts model months to `Dates.Date`s (first day of the month), the time axis expected
by PlanetTools. Only year and month of the dates are used by PlanetTools, so the
Gregorian calendar of `Date` does not affect results for other planets.

# Example

```julia
julia> yearmonth_dates([YearMonth(3, 2)])
1-element Vector{Dates.Date}:
 0003-02-01
```
"""
yearmonth_dates(months::AbstractVector{YearMonth}) = [Date(m.year, m.month) for m in months]

"""
    annual_means(x, sm, months, month_weights; threshold=80.0) -> PeriodMeans

Annual means of all complete calendar years of the monthly record, for the gridded field
`x` and its zonal and global means `sm`, computed with
`PlanetTools.PlanetStatistics.calc_annual_mean_ts!`.

Definitions (PlanetTools): month `t` gets the weight `month_weights[t]` times the
availability of its input (1 for the gridded field; `sm.zm_w`, `sm.gm_w` for the zonal
and global means, i.e. months with partly missing spatial coverage get less weight); the
mean is `NaN` if fewer than `threshold` percent of the 12 months are valid; the
availability of the annual mean is the fraction of valid months.

# Arguments

  - `x::AbstractArray{<:Real}`: monthly field `(Lon, Lat, [Lev], Time)`.
  - `sm::SpatialMeans`: its zonal and global means from [`planet_means`](@ref).
  - `months::AbstractVector{YearMonth}`: labels of the `Time` dimension (consecutive).
  - `month_weights`: weights of the months (length `Time`), e.g.
    [`month_weights`](@ref)`(calendar, months)`, or `nothing` for equal weights.

# Keyword arguments

  - `threshold::Real = 80.0`: minimum percentage of valid months.

# Returns

[`PeriodMeans`](@ref) with one entry per year of [`annual_periods`](@ref)`(months)`; an
empty last dimension if there is no complete year.

# Example

```julia
w = month_weights(MARS_CALENDAR, months)
an = annual_means(x, planet_means(x, area), months, w)
```
"""
function annual_means(
    x::AbstractArray{<:Real},
    sm::SpatialMeans,
    months::AbstractVector{YearMonth},
    month_weights;
    threshold::Real=80.0,
)
    years = annual_periods(months)
    out = alloc_period_means(x, sm, length(years))
    isempty(years) && return out
    full_years = Date(first(years).first.year):Year(1):Date(last(years).first.year)
    calc_annual_mean_ts!(
        lift(x),
        lift(sm.zm),
        lift(sm.gm),
        lift(out.ij),
        lift(out.zm),
        lift(out.gm),
        lift(out.ij_w),
        lift(out.zm_w),
        lift(out.gm_w),
        yearmonth_dates(months),
        [length(months)],
        [full_years];
        zm_w=lift(sm.zm_w),
        gm_w=lift(sm.gm_w),
        month_weights=month_weights,
        threshold=threshold,
    )
    return out
end

"""
    seasonal_means(x, sm, months, season, month_weights; threshold=80.0) -> PeriodMeans

Means over all complete occurrences of `season` in the monthly record, for the gridded
field `x` and its zonal and global means `sm`, computed with
`PlanetTools.PlanetStatistics.calc_seasonal_mean_ts!`. Weighting, availability, and
threshold are defined as for [`annual_means`](@ref), with the months of the season in
place of the 12 months of a year.

# Arguments

  - `x`, `sm`, `months`: as for [`annual_means`](@ref).
  - `season::SeasonSpec`: the season.
  - `month_weights`: weights of the months, or `nothing` for equal weights.

# Keyword arguments

  - `threshold::Real = 80.0`: minimum percentage of valid months of a season.

# Returns

[`PeriodMeans`](@ref) with one entry per element of
[`seasonal_periods`](@ref)`(months, season)`.

# Example

```julia
se = seasonal_means(x, sm, months, SeasonSpec("FMA", 2, 3), w)
```
"""
function seasonal_means(
    x::AbstractArray{<:Real},
    sm::SpatialMeans,
    months::AbstractVector{YearMonth},
    season::SeasonSpec,
    month_weights;
    threshold::Real=80.0,
)
    periods = seasonal_periods(months, season)
    np = length(periods)
    out = alloc_period_means(x, sm, np)
    np == 0 && return out
    ind = reshape([first(p.slots) for p in periods], 1, 1, np)   # (Sim, Season, Year)
    l3(a) = reshape(a, 1, 1, 1, size(a)...)                       # (Sim, Var, Season, …)
    calc_seasonal_mean_ts!(
        lift(x),
        lift(sm.zm),
        lift(sm.gm),
        l3(out.ij),
        l3(out.zm),
        l3(out.gm),
        l3(out.ij_w),
        l3(out.zm_w),
        l3(out.gm_w),
        ind,
        yearmonth_dates(months);
        zm_w=lift(sm.zm_w),
        gm_w=lift(sm.gm_w),
        month_weights=month_weights,
        length_of_season=season.length,
        threshold=threshold,
    )
    return out
end

"""
    alloc_period_means(x, sm, np) -> PeriodMeans

`NaN32`-initialized [`PeriodMeans`](@ref) for `np` periods, with the spatial dimensions
of `x` and `sm`.

# Example

```julia
pm = alloc_period_means(x, sm, 20)
```
"""
function alloc_period_means(x, sm::SpatialMeans, np::Integer)
    sp(a) = size(a)[1:(end - 1)]
    return PeriodMeans(
        nanarray(sp(x)..., np),
        nanarray(sp(x)..., np),
        nanarray(sp(sm.zm)..., np),
        nanarray(sp(sm.zm)..., np),
        nanarray(sp(sm.gm)..., np),
        nanarray(sp(sm.gm)..., np),
    )
end

"""
    multiyear_stats(x, x_w, slots; threshold=80.0, corrected=true) -> (mean, sd, avail)

Multi-year mean and standard deviation over the contiguous samples `slots` (last
dimension of `x`, e.g. the annual or seasonal means of the years in the statistics
window), computed with `PlanetTools.PlanetStatistics.calc_multiyear_stats!`
(`Val(:annual)`).

Definitions (PlanetTools): with the availabilities `x_w` as weights `w` of the valid
samples, `V₁ = Σw`, `V₂ = Σw²`, the mean is `μ = Σwx / V₁` and the standard deviation

  - `corrected=true` (default): `σ = √(Σw(x − μ)² / (V₁ − V₂/V₁))`, the sample standard
    deviation for reliability weights, `√(Σ(x − μ)² / (n − 1))` for equal weights; `NaN`
    with fewer than two valid samples;
  - `corrected=false`: the population standard deviation `√(Σw(x − μ)² / V₁)`.

Results are `NaN` if fewer than `threshold` percent of the samples are valid; `avail` is
the fraction of valid samples.

# Arguments

  - `x::AbstractArray{<:Real}`: samples along the last dimension.
  - `x_w`: availability weights of `x` (same size) or `nothing` (equal weights).
  - `slots::AbstractUnitRange{<:Integer}`: indices of the samples to use.

# Keyword arguments

  - `threshold::Real = 80.0`: minimum percentage of valid samples.
  - `corrected::Bool = true`: sample (`true`) or population (`false`) standard deviation.

# Example

```julia
julia> me, sd, av = multiyear_stats(reshape(Float32[1, 2, 3, NaN], 1, 4), nothing, 1:4;
                                    threshold=50);

julia> me[1], sd[1], av[1]
(2.0f0, 1.0f0, 0.75f0)
```
"""
function multiyear_stats(
    x::AbstractArray{<:Real},
    x_w,
    slots::AbstractUnitRange{<:Integer};
    threshold::Real=80.0,
    corrected::Bool=true,
)
    sp = size(x)[1:(end - 1)]
    me, sd, av = nanarray(sp...), nanarray(sp...), nanarray(sp...)
    isempty(slots) && return me, sd, av
    calc_multiyear_stats!(
        lift(me),
        lift(sd),
        lift(av),
        nothing,
        nothing,
        lift(x),
        [first(slots) last(slots)],
        Val(:annual);
        data_w=lift(x_w),
        threshold=threshold,
        corrected=corrected,
    )
    return me, sd, av
end

"""
    monthly_stats(x, x_w, slots12; threshold=80.0, corrected=true) -> (mean, sd, avail)

Multi-year mean and standard deviation of each calendar month of the annual cycle,
computed with `PlanetTools.PlanetStatistics.calc_multiyear_stats!` (`Val(:monthly)`),
with the same definitions as [`multiyear_stats`](@ref) applied month by month.

# Arguments

  - `x::AbstractArray{<:Real}`: monthly field, time as last dimension.
  - `x_w`: availability weights of `x` (e.g. `sm.zm_w` for zonal means), or `nothing`.
  - `slots12::AbstractVector{<:AbstractVector{<:Integer}}`: 12 index vectors, see
    [`climatology_slots`](@ref).

# Keyword arguments

  - `threshold::Real = 80.0`: minimum percentage of valid samples per month.
  - `corrected::Bool = true`: sample (`true`) or population (`false`) standard deviation.

# Returns

Arrays with the leading dimensions of `x` and 12 as last dimension.

# Example

```julia
me, sd, av = monthly_stats(sm.gm, sm.gm_w, climatology_slots(months, window))
```
"""
function monthly_stats(
    x::AbstractArray{<:Real},
    x_w,
    slots12::AbstractVector{<:AbstractVector{<:Integer}};
    threshold::Real=80.0,
    corrected::Bool=true,
)
    length(slots12) == 12 || throw(ArgumentError("slots12 must have 12 elements"))
    sp = size(x)[1:(end - 1)]
    me, sd, av = nanarray(sp..., 12), nanarray(sp..., 12), nanarray(sp..., 12)
    idx = reshape([collect(Int, s) for s in slots12], 1, 12)
    calc_multiyear_stats!(
        lift(me),
        lift(sd),
        lift(av),
        nothing,
        nothing,
        lift(x),
        idx,
        Val(:monthly);
        data_w=lift(x_w),
        threshold=threshold,
        corrected=corrected,
    )
    return me, sd, av
end

"""
    contiguous_range(idx::AbstractVector{<:Integer}) -> UnitRange{Int}

Converts a vector of consecutive indices into a `UnitRange` (as required by
`calc_multiyear_stats!` with `Val(:annual)`); throws an `ArgumentError` if the indices are
not consecutive. An empty vector gives an empty range.

# Example

```julia
julia> contiguous_range([3, 4, 5])
3:5
```
"""
function contiguous_range(idx::AbstractVector{<:Integer})
    isempty(idx) && return 1:0
    r = first(idx):last(idx)
    collect(r) == idx || throw(ArgumentError("indices are not contiguous: $idx"))
    return r
end
