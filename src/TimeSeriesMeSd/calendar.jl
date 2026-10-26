# Calendar handling for planetary (Earth, Mars, exoplanet) model calendars.
#
# ModelE/ROCKE-3D organize output into 12 "months" per model year. On planets other than
# Earth the months have different lengths (in local solar days). Only the year/month
# label and the month lengths are needed for the statistics, so a light-weight
# `YearMonth` type is used instead of `Dates.Date`, which would impose the Gregorian
# calendar and does not allow, e.g., a 56-sol month.

const DEFAULT_MONTH_NAMES = (
    "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"
)

"""
    YearMonth(year::Integer, month::Integer)

Label of one model month: model year `year` (any integer, typically ≥ 1 for ModelE and
ROCKE-3D runs) and month `month` (1 … 12). `YearMonth`s are totally ordered (earlier <
later), can be shifted by an integer number of months with `+`/`-`, and the number of
months between two labels is obtained with `months_between`.

The type deliberately carries no information on the length of a month; that is a
property of the [`PlanetCalendar`](@ref).

# Arguments

  - `year::Integer`: model year.
  - `month::Integer`: month of the year, `1 ≤ month ≤ 12`; anything else throws an
    `ArgumentError`.

# Example

```julia
julia> a = YearMonth(3, 11)
YearMonth(3, 11)

julia> a + 3
YearMonth(4, 2)

julia> months_between(YearMonth(1, 1), a)
34

julia> YearMonth(1, 12) < YearMonth(2, 1)
true
```
"""
struct YearMonth
    year::Int
    month::Int
    function YearMonth(year::Integer, month::Integer)
        1 <= month <= 12 ||
            throw(ArgumentError("month must be between 1 and 12, got $month"))
        return new(Int(year), Int(month))
    end
end

"""
    absmonth(ym::YearMonth) -> Int

Absolute month counter `12 * year + month - 1`, used for ordering and arithmetic of
[`YearMonth`](@ref) values. Consecutive months have consecutive counters.

# Example

```julia
julia> absmonth(YearMonth(0, 1)), absmonth(YearMonth(1, 1))
(0, 12)
```
"""
absmonth(ym::YearMonth) = 12 * ym.year + ym.month - 1

"""
    fromabsmonth(n::Integer) -> YearMonth

Inverse of [`absmonth`](@ref): converts an absolute month counter back into a
[`YearMonth`](@ref). Works for negative counters (years before year 0) as well.

# Example

```julia
julia> fromabsmonth(35)
YearMonth(2, 12)
```
"""
fromabsmonth(n::Integer) = YearMonth(fld(n, 12), mod(n, 12) + 1)

Base.isless(a::YearMonth, b::YearMonth) = absmonth(a) < absmonth(b)
Base.:+(a::YearMonth, n::Integer) = fromabsmonth(absmonth(a) + n)
Base.:-(a::YearMonth, n::Integer) = fromabsmonth(absmonth(a) - n)
Base.show(io::IO, a::YearMonth) = print(io, "YearMonth(", a.year, ", ", a.month, ")")

"""
    months_between(a::YearMonth, b::YearMonth) -> Int

Number of months from `a` to `b` (`b - a`); negative if `b` precedes `a`.

# Example

```julia
julia> months_between(YearMonth(1, 1), YearMonth(20, 12))
239
```
"""
months_between(a::YearMonth, b::YearMonth) = absmonth(b) - absmonth(a)

"""
    yearmonth_range(first::YearMonth, last::YearMonth) -> Vector{YearMonth}

All months from `first` to `last`, both inclusive, in chronological order. Returns an
empty vector if `last < first`.

# Example

```julia
julia> yearmonth_range(YearMonth(1, 11), YearMonth(2, 2))
4-element Vector{YearMonth}:
 YearMonth(1, 11)
 YearMonth(1, 12)
 YearMonth(2, 1)
 YearMonth(2, 2)
```
"""
function yearmonth_range(first::YearMonth, last::YearMonth)
    return [first + k for k in 0:months_between(first, last)]
end

"""
    yyyy(year::Integer) -> String

Four-digit, zero-padded representation of a model year as used in ModelE file names and
in the output file names of this package (`3 → "0003"`).

# Example

```julia
julia> yyyy(3)
"0003"
```
"""
yyyy(year::Integer) = lpad(string(year), 4, '0')

"""
    PlanetCalendar(name, days_in_month; day_length="", month_names=DEFAULT_MONTH_NAMES)

Model calendar of a planet: the lengths of the 12 model months and descriptive metadata.

The month lengths are used as weights when monthly means are averaged into seasonal and
annual means: month `m` gets the weight `days_in_month[m] / days_in_year(calendar)`. The
unit of `days_in_month` is irrelevant for the weights (only ratios enter), but by
convention it is the number of local solar days (sols on Mars).

# Arguments

  - `name::AbstractString`: planet or calendar name, written to the output metadata.
  - `days_in_month`: 12 month lengths (any real numbers > 0), in local days.

# Keyword arguments

  - `day_length::AbstractString = ""`: length of one local day with unit, e.g.
    `"88775.245 seconds"` for Mars; written as attribute `day_length` of the time
    coordinate (omitted when empty).
  - `month_names`: 12 three-letter month labels used in output file names.
  - `cf_calendar::AbstractString = lowercase(name)`: value of the CF `calendar` attribute
    of the time coordinate. A name that is not a CF standard calendar (e.g. `"mars"`)
    declares a user-defined calendar, whose month lengths are written as attribute
    `month_lengths`; CF readers then do not convert the time values into Earth dates.

# Predefined calendars

  - [`MARS_CALENDAR`](@ref): ROCKE-3D Mars calendar, months of 61, 66, 66, 65, 60, 54, 50,
    47, 46, 47, 51, 56 sols (669 sols per year).
  - [`EARTH_NOLEAP_CALENDAR`](@ref): ModelE 365-day calendar without leap years.

# Example

```julia
julia> cal = PlanetCalendar("TOI-700d", fill(10.0, 12); day_length="864000 seconds");

julia> days_in_year(cal)
120.0
```
"""
struct PlanetCalendar
    name::String
    days_in_month::NTuple{12,Float64}
    day_length::String
    month_names::NTuple{12,String}
    cf_calendar::String
end

function PlanetCalendar(
    name::AbstractString,
    days_in_month;
    day_length::AbstractString="",
    month_names=DEFAULT_MONTH_NAMES,
    cf_calendar::AbstractString=lowercase(name),
)
    length(days_in_month) == 12 ||
        throw(ArgumentError("days_in_month must have 12 elements"))
    all(>(0), days_in_month) || throw(ArgumentError("month lengths must be positive"))
    length(month_names) == 12 || throw(ArgumentError("month_names must have 12 elements"))
    return PlanetCalendar(
        String(name),
        NTuple{12,Float64}(Float64.(collect(days_in_month))),
        String(day_length),
        NTuple{12,String}(String.(collect(month_names))),
        String(cf_calendar),
    )
end

"""
    MARS_CALENDAR

[`PlanetCalendar`](@ref) of the ROCKE-3D Mars configuration: 12 months of 61, 66, 66, 65,
60, 54, 50, 47, 46, 47, 51, 56 sols (669 sols per Mars year), and a sol of
88775.245 seconds.

# Example

```julia
julia> month_weights(MARS_CALENDAR, [YearMonth(1, 2)])
1-element Vector{Float64}:
 0.09865470852017937
```
"""
const MARS_CALENDAR = PlanetCalendar(
    "Mars", (61, 66, 66, 65, 60, 54, 50, 47, 46, 47, 51, 56); day_length="88775.245 seconds"
)

"""
    EARTH_NOLEAP_CALENDAR

[`PlanetCalendar`](@ref) of ModelE for Earth: 365-day year without leap years (months of
31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 days) and a day of 86400 seconds.

# Example

```julia
julia> days_in_year(EARTH_NOLEAP_CALENDAR)
365.0
```
"""
const EARTH_NOLEAP_CALENDAR = PlanetCalendar(
    "Earth",
    (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);
    day_length="86400 seconds",
    cf_calendar="noleap",
)

"""
    days_in_year(cal::PlanetCalendar) -> Float64

Sum of the 12 month lengths of the calendar.

# Example

```julia
julia> days_in_year(MARS_CALENDAR)
669.0
```
"""
days_in_year(cal::PlanetCalendar) = sum(cal.days_in_month)

"""
    month_weights(cal::PlanetCalendar, months::AbstractVector{YearMonth}) -> Vector{Float64}

Weight of every month in `months` for temporal averaging: the month length divided by
the length of the year. For a complete year the weights add up to one.

# Example

```julia
julia> w = month_weights(MARS_CALENDAR, yearmonth_range(YearMonth(1, 1), YearMonth(1, 12)));

julia> sum(w) ≈ 1
true
```
"""
function month_weights(cal::PlanetCalendar, months::AbstractVector{YearMonth})
    diy = days_in_year(cal)
    return [cal.days_in_month[ym.month] / diy for ym in months]
end

"""
    monthname(cal::PlanetCalendar, ym::YearMonth) -> String

Three-letter label of the month of `ym` in calendar `cal` (e.g. `"JAN"`).

# Example

```julia
julia> monthname(MARS_CALENDAR, YearMonth(4, 12))
"DEC"
```
"""
monthname(cal::PlanetCalendar, ym::YearMonth) = cal.month_names[ym.month]

"""
    SeasonSpec(name::AbstractString, first_month::Integer, length::Integer)

A season of `length` consecutive months starting in month `first_month`. Seasons may
wrap into the following year (e.g. `SeasonSpec("DJF", 12, 3)`); the season is then
labeled with the year of its *first* month (the year of the December). `name` is used
in output file names, e.g. `FMA_0001-0020.<name><run>.ij.sts.nc`.

# Arguments

  - `name`: label of the season (e.g. `"JJA"`).
  - `first_month`: first month, `1 ≤ first_month ≤ 12`.
  - `length`: number of months, `1 ≤ length ≤ 12`.

# Predefined season sets

  - [`EARTH_SEASONS`](@ref): DJF, MAM, JJA, SON.
  - [`MARS_SEASONS`](@ref): FMA, MJJ, ASO, NDJ.

# Example

```julia
julia> s = SeasonSpec("NDJ", 11, 3)
SeasonSpec("NDJ", 11, 3)

julia> season_months(s)
3-element Vector{Int64}:
 11
 12
  1
```
"""
struct SeasonSpec
    name::String
    first_month::Int
    length::Int
    function SeasonSpec(name::AbstractString, first_month::Integer, length::Integer)
        1 <= first_month <= 12 ||
            throw(ArgumentError("first_month must be in 1:12, got $first_month"))
        1 <= length <= 12 || throw(ArgumentError("length must be in 1:12, got $length"))
        return new(String(name), Int(first_month), Int(length))
    end
end

function Base.show(io::IO, s::SeasonSpec)
    return print(io, "SeasonSpec(\"", s.name, "\", ", s.first_month, ", ", s.length, ")")
end

"""
    season_months(s::SeasonSpec) -> Vector{Int}

Calendar months (1 … 12) belonging to season `s`, in chronological order.

# Example

```julia
julia> season_months(SeasonSpec("DJF", 12, 3))
3-element Vector{Int64}:
 12
  1
  2
```
"""
season_months(s::SeasonSpec) = [mod1(s.first_month + k, 12) for k in 0:(s.length - 1)]

"""
    EARTH_SEASONS

Standard meteorological seasons for Earth: `DJF`, `MAM`, `JJA`, `SON`.

# Example

```julia
julia> [s.name for s in EARTH_SEASONS]
4-element Vector{String}:
 "DJF"
 "MAM"
 "JJA"
 "SON"
```
"""
const EARTH_SEASONS = [
    SeasonSpec("DJF", 12, 3),
    SeasonSpec("MAM", 3, 3),
    SeasonSpec("JJA", 6, 3),
    SeasonSpec("SON", 9, 3),
]

"""
    MARS_SEASONS

Season set for ROCKE-3D Mars: `FMA`, `MJJ`, `ASO`, `NDJ` (three model months each,
starting in February, May, August, November).

# Example

```julia
julia> [s.first_month for s in MARS_SEASONS]
4-element Vector{Int64}:
  2
  5
  8
 11
```
"""
const MARS_SEASONS = [
    SeasonSpec("FMA", 2, 3),
    SeasonSpec("MJJ", 5, 3),
    SeasonSpec("ASO", 8, 3),
    SeasonSpec("NDJ", 11, 3),
]

"""
    TimeAxis(reference::YearMonth, units::AbstractString)
    TimeAxis(; reference=YearMonth(0, 12), units="months since 0000-12")

Definition of the numerical time coordinate written to the time-series output files.

The coordinate value of a period (month, season, or year) is the number of months from
`reference` to the *first* month of that period, `months_between(reference, first)`.
One common unit (months) is therefore used for monthly, seasonal, and annual series.
The string `units` is written verbatim as the `units` attribute and must describe
`reference` in the convention of the user.

With the default `reference = YearMonth(0, 12)`, the first month of model year 1
(`YearMonth(1, 1)`) has the coordinate value 1, so a monthly series starting in January
of year 1 has `time = 1, 2, …, n`.

In addition, every time-series file contains the integer variables `year(time)` and
`month(time)` with the first model month of each period, which are unambiguous
irrespective of the units string.

# Example

```julia
julia> ax = TimeAxis(; reference=YearMonth(0, 12), units="months since 0023-12-56");

julia> time_value(ax, YearMonth(2, 1))
13
```
"""
struct TimeAxis
    reference::YearMonth
    units::String
end

function TimeAxis(;
    reference::YearMonth=YearMonth(0, 12), units::AbstractString="months since 0000-12"
)
    return TimeAxis(reference, String(units))
end

"""
    time_value(ax::TimeAxis, first_month::YearMonth) -> Int

Numerical time coordinate of a period starting in `first_month`, see [`TimeAxis`](@ref).

# Example

```julia
julia> time_value(TimeAxis(), YearMonth(1, 1))
1
```
"""
time_value(ax::TimeAxis, first_month::YearMonth) = months_between(ax.reference, first_month)

## Period definitions -------------------------------------------------------------------

"""
    Period

An averaging period in a monthly record: `first` is the first month (label), `slots` the
indices of all months of the period in the monthly record (consecutive, ascending).

# Example

```julia
julia> p = Period(YearMonth(1, 1), 1:12);

julia> length(p.slots)
12
```
"""
struct Period
    first::YearMonth
    slots::UnitRange{Int}
end

"""
    annual_periods(months::AbstractVector{YearMonth}) -> Vector{Period}

All complete calendar years (January … December) contained in the monthly record
`months`, which must be consecutive months in chronological order. Partial years at the
beginning or end of the record are omitted. Each [`Period`](@ref) holds the label of the
January and the 12 record indices of the year.

# Example

```julia
julia> months = yearmonth_range(YearMonth(1, 7), YearMonth(3, 12));

julia> [(p.first.year, p.slots) for p in annual_periods(months)]
2-element Vector{Tuple{Int64, UnitRange{Int64}}}:
 (2, 7:18)
 (3, 19:30)
```
"""
function annual_periods(months::AbstractVector{YearMonth})
    return seasonal_periods(months, SeasonSpec("ANN", 1, 12))
end

"""
    seasonal_periods(months::AbstractVector{YearMonth}, season::SeasonSpec) -> Vector{Period}

All complete occurrences of `season` in the monthly record `months` (consecutive months
in chronological order). An occurrence is complete if all `season.length` months are
inside the record; seasons wrapping into the next year (e.g. NDJ) are therefore only
included if the January of the following year is part of the record.

# Example

```julia
julia> months = yearmonth_range(YearMonth(1, 1), YearMonth(2, 12));

julia> [(p.first, p.slots) for p in seasonal_periods(months, SeasonSpec("NDJ", 11, 3))]
1-element Vector{Tuple{YearMonth, UnitRange{Int64}}}:
 (YearMonth(1, 11), 11:13)
```
"""
function seasonal_periods(months::AbstractVector{YearMonth}, season::SeasonSpec)
    check_consecutive(months)
    periods = Period[]
    n = length(months)
    for t in 1:(n - season.length + 1)
        if months[t].month == season.first_month
            push!(periods, Period(months[t], t:(t + season.length - 1)))
        end
    end
    return periods
end

"""
    period_in_window(p::Period, months, window::Tuple{YearMonth,YearMonth}) -> Bool

`true` if all months of period `p` (indices into the monthly record `months`) lie inside
the closed interval `window = (first, last)`.

# Example

```julia
julia> months = yearmonth_range(YearMonth(1, 1), YearMonth(3, 12));

julia> ps = annual_periods(months);

julia> [period_in_window(p, months, (YearMonth(2, 1), YearMonth(3, 12))) for p in ps]
3-element Vector{Bool}:
 0
 1
 1
```
"""
function period_in_window(
    p::Period, months::AbstractVector{YearMonth}, window::Tuple{YearMonth,YearMonth}
)
    return window[1] <= months[first(p.slots)] && months[last(p.slots)] <= window[2]
end

"""
    climatology_slots(months, window) -> Vector{Vector{Int}}

Record indices of every calendar month (1 … 12) to be used for the multi-year statistics
of the annual cycle (`mom`, `som` products). Only *complete calendar years* that lie
entirely inside `window` and inside the record are used, so that each calendar month is
sampled the same number of times.

Returns a vector of 12 index vectors; element `m` contains the indices of month `m`.

# Example

```julia
julia> months = yearmonth_range(YearMonth(1, 1), YearMonth(4, 6));

julia> slots = climatology_slots(months, (YearMonth(2, 3), YearMonth(4, 6)));

julia> slots[1]   # only year 3 is a full year inside the window
1-element Vector{Int64}:
 25
```
"""
function climatology_slots(
    months::AbstractVector{YearMonth}, window::Tuple{YearMonth,YearMonth}
)
    slots = [Int[] for _ in 1:12]
    for p in annual_periods(months)
        period_in_window(p, months, window) || continue
        for t in p.slots
            push!(slots[months[t].month], t)
        end
    end
    return slots
end

"""
    check_consecutive(months::AbstractVector{YearMonth})

Throws an `ArgumentError` unless `months` consists of consecutive months in chronological
order (each element one month after its predecessor).

# Example

```julia
julia> check_consecutive(yearmonth_range(YearMonth(1, 1), YearMonth(1, 6)))  # no error
```
"""
function check_consecutive(months::AbstractVector{YearMonth})
    for t in 2:length(months)
        months_between(months[t - 1], months[t]) == 1 ||
            throw(ArgumentError("monthly record is not consecutive at index $t"))
    end
    return nothing
end
