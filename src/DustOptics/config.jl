# Configuration of the optical post-processing.

"""
    OPTICAL_QUANTITIES

Optical quantities that can be computed, with the prefix of their variable names and the
first word of their long names: `:sct` (scattering), `:abs` (absorption), and `:ext`
(extinction = scattering + absorption).

# Example

```julia
julia> OPTICAL_QUANTITIES[:ext]
("ext", "Extinction")
```
"""
const OPTICAL_QUANTITIES = OrderedDict(
    :sct => ("sct", "Scattering"),
    :abs => ("abs", "Absorption"),
    :ext => ("ext", "Extinction"),
)

"""
    DustOpticsConfig(; stats, bins, input_dir, input_name, kwargs...)

All settings of the optical post-processing of dust, [`run_dust_optics`](@ref).

The input are the monthly time series of the dust tracer mass mixing ratios and of the
air mass of the layers written by `PlanetPostStats.TimeSeriesMeSd`
(`<prefix>.<input_name><run>.<input_tag>.mts.nc` in
`<input_root>/<run>/<input_dir>/MonthlyTs/`). From these, for every run and month, the
dust masses of the size `bins` in the model layers are formed, and from them

  - the effective radius of the dust in the layers (`effRad3d_<class>`), of the column
    dust load (`effRadColMass_<class>`), and its mass-weighted column average
    (`effRadColAve_<class>`), where `<class>` is the largest size class of
    [`optics_classes`](@ref) (e.g. `PM32Dust`), or `Dust` for other bins;
  - for every [`SpectralSet`](@ref), diagnostic wavelength, and quantity of
    `quantities`, the optical depth of every bin and every size class in the layers
    (`<q>3d_<λ>_<bin>`, e.g. `ext3d_550nm_Clay1`) and of the column (`<q>2d_<λ>_<bin>`,
    e.g. `ext2d_550nm_PM32Dust`).

All variables are then processed like the variables of `TimeSeriesMeSd` (time series,
multi-year means and standard deviations; gridded, zonal, and global means), with the
settings of `stats`.

# Required keyword arguments

  - `stats::TimeSeriesConfig`: runs, periods, calendar, seasons, time axis, vertical
    coordinate, optional pressure levels, threshold, products, and output settings
    (`simpath`, `output_dir`, `output_name`, `global_attributes`, …). Its `groups` and
    `derived` are not used and should be empty.
  - `bins::Vector{DustBin}`: size bins of the optical calculations, in the order of
    accumulation of the size classes.
  - `input_dir::String`: subdirectory of a run with the `TimeSeriesMeSd` output of the
    dust tracers (its `output_dir`, e.g. `"Tracermass3d"`).
  - `input_name::String`: tag of these files (its `output_name`, e.g. `"tracermass3d"`).

# Optional keyword arguments

  - `spectral::Vector{SpectralSet} = SpectralSet[]`: spectral ranges, wavelengths, and
    lookup tables; empty: only effective radii.
  - `fractions = nothing`: fractions of the source tracer mass in the bins, one value per
    bin: a vector (all runs), a `Dict` from run identifier to vector, or a function
    `run -> vector`; `nothing`: all fractions 1. See [`bin_fractions`](@ref). The
    fractions of the bins of one source should sum to 1.
  - `input_root = nothing`: directory with one subdirectory per run for the input;
    `nothing`: `stats.simpath`.
  - `input_tag::String = "ijl"`: representation tag of the input files (gridded, model
    levels).
  - `input_prefix = nothing`: period prefix of the input file names (e.g.
    `"JAN0001-DEC0020"`); `nothing`: the single file matching
    `*.<input_name><run>.<input_tag>.mts.nc`.
  - `airmass::String = "airmass"`: input variable with the air mass of the layers per
    unit area (kg m⁻²).
  - `quantities::Vector{Symbol} = [:sct, :abs, :ext]`: optical quantities, see
    [`OPTICAL_QUANTITIES`](@ref).
  - `layers::Bool = true`: output the optical depths of the layers (3D).
  - `columns::Bool = true`: output the column optical depths (2D).
  - `bin_variables::Bool = true`: output the optical depths of the individual bins.
  - `class_variables::Bool = true`: output the optical depths of the size classes.
  - `extrapolation::Symbol = :linear`: extrapolation of the lookup tables, see
    [`interpolate_coefficients`](@ref).
  - `effective_radius_to_pressure::Bool = true`: interpolate the effective radius of the
    layers to pressure levels when `stats.pressure_levels` is set. The optical depths of
    the layers are layer-integrated quantities and are never interpolated.

With `stats.pressure_levels`, the surface pressure (`stats.pressure_levels.
surface_pressure`, e.g. `prsurf`) must be contained in the input files; it is read but
not written.

# Example

```julia
stats = TimeSeriesConfig(; simpath="/path/to/simulations", runs=["MarsDustRun1"],
    sim_period=(YearMonth(1, 1), YearMonth(20, 12)),
    stat_period=(YearMonth(3, 1), YearMonth(20, 12)), calendar=MARS_CALENDAR,
    seasons=MARS_SEASONS, output_dir="Traceroptics3d", output_name="traceroptics3d")
bins = [DustBin("Clay1"; source="Clay", radius=0.132e-6, density=2.5e3), …]
cfg = DustOpticsConfig(; stats=stats, bins=bins, spectral=[sw, lw],
    fractions=bin_fractions(bins, "Clay" => [0.009, 0.081, 0.234, 0.676]),
    input_dir="Tracermass3d", input_name="tracermass3d")
```
"""
Base.@kwdef struct DustOpticsConfig
    stats::TimeSeriesConfig
    bins::Vector{DustBin}
    input_dir::String
    input_name::String
    spectral::Vector{SpectralSet} = SpectralSet[]
    fractions::Any = nothing
    input_root::Union{Nothing,String} = nothing
    input_tag::String = "ijl"
    input_prefix::Union{Nothing,String} = nothing
    airmass::String = "airmass"
    quantities::Vector{Symbol} = [:sct, :abs, :ext]
    layers::Bool = true
    columns::Bool = true
    bin_variables::Bool = true
    class_variables::Bool = true
    extrapolation::Symbol = :linear
    effective_radius_to_pressure::Bool = true
end

"""
    run_fractions(cfg::DustOpticsConfig, run) -> Vector{Float64}

Fractions of the source tracer mass in the bins for run `run`, from `cfg.fractions`
(vector, `Dict` by run, function of the run, or `nothing` for all ones). Throws an
`ArgumentError` if the number of fractions differs from the number of bins, a fraction
is negative or not finite, or a `Dict` has no entry for the run; warns if the fractions
of the bins of one source do not sum to 1 (within 10⁻³).

# Example

```julia
julia> run_fractions(cfg, "MarsDustRun1")
8-element Vector{Float64}: …
```
"""
function run_fractions(cfg::DustOpticsConfig, run::AbstractString)
    fr = cfg.fractions
    f = if fr === nothing
        ones(length(cfg.bins))
    elseif fr isa AbstractVector
        fr
    elseif fr isa AbstractDict
        haskey(fr, run) || throw(ArgumentError("fractions: no entry for run $run"))
        fr[run]
    else
        fr(String(run))
    end
    f = Float64.(collect(f))
    length(f) == length(cfg.bins) || throw(
        ArgumentError("run $run: $(length(f)) fractions for $(length(cfg.bins)) bins")
    )
    all(x -> isfinite(x) && x >= 0, f) ||
        throw(ArgumentError("run $run: fractions must be non-negative and finite"))
    for s in bin_sources(cfg.bins)
        total = sum(f[i] for (i, b) in enumerate(cfg.bins) if b.source == s)
        abs(total - 1) <= 1e-3 ||
            @warn "Run $run: the fractions of the bins of source $s sum to $total"
    end
    return f
end

"""
    effective_radius_suffix(bins) -> String

Suffix of the names of the effective-radius variables: the largest size class of
[`optics_classes`](@ref)`(bins)` (e.g. `"PM32Dust"`), or `"Dust"` if there is none.

# Example

```julia
julia> effective_radius_suffix([DustBin("x"; radius=1e-6, density=2500)])
"Dust"
```
"""
function effective_radius_suffix(bins::AbstractVector{DustBin})
    classes = optics_classes(bins)
    return isempty(classes) ? "Dust" : last(classes).first
end

"""
    optics_variable_names(cfg::DustOpticsConfig) -> Vector{String}

Names of all output variables of the optical post-processing, in output order: the
effective radii, then, for every spectral set, quantity, wavelength, and for layers and
columns, the bins followed by the size classes.

# Example

```julia
julia> optics_variable_names(cfg)[1:3]
3-element Vector{String}:
 "effRad3d_PM32Dust"
 "effRadColMass_PM32Dust"
 "effRadColAve_PM32Dust"
```
"""
function optics_variable_names(cfg::DustOpticsConfig)
    return [spec.name for spec in optics_variable_specs(cfg, nothing)]
end

"""
    OpticsVariableSpec

Description of one output variable of the optical post-processing (internal).

# Fields

  - `name::String`, `long_name::String`, `units::String`.
  - `kind::Symbol`: `:reff3d`, `:reffcolmass`, `:reffcolave`, or `:tau`.
  - `is3d::Bool`: layer (3D) or column (2D) variable.
  - `coefficients::Vector{Float64}`: for `:tau`, coefficient of every bin
    (m² kg⁻¹; zero for bins not included).

# Example

```julia
specs = optics_variable_specs(cfg, coefficients)
first(specs).name                          # "effRad3d_PM32Dust"
```
"""
struct OpticsVariableSpec
    name::String
    long_name::String
    units::String
    kind::Symbol
    is3d::Bool
    coefficients::Vector{Float64}
end

"""
    optics_variable_specs(cfg::DustOpticsConfig, coefficients) -> Vector{OpticsVariableSpec}

Descriptions of all output variables (see [`optics_variable_names`](@ref)).
`coefficients` are the [`OpticalCoefficients`](@ref) of the spectral sets, in the order
of `cfg.spectral`; with `nothing`, only names and attributes are formed (zero
coefficients), e.g. for validation.

# Example

```julia
coefs = [interpolate_coefficients(s, cfg.bins) for s in cfg.spectral]
specs = optics_variable_specs(cfg, coefs)
```
"""
function optics_variable_specs(cfg::DustOpticsConfig, coefficients)
    bins = cfg.bins
    nb = length(bins)
    sfx = effective_radius_suffix(bins)
    specs = OpticsVariableSpec[
        OpticsVariableSpec(
            "effRad3d_" * sfx, sfx * " Effective Radius", "m", :reff3d, true, Float64[]
        ),
        OpticsVariableSpec(
            "effRadColMass_" * sfx,
            sfx * " Effective Radius Column Mass",
            "m",
            :reffcolmass,
            false,
            Float64[],
        ),
        OpticsVariableSpec(
            "effRadColAve_" * sfx,
            sfx * " Effective Radius Column Average",
            "m",
            :reffcolave,
            false,
            Float64[],
        ),
    ]
    classes = optics_classes(bins)
    for (s, set) in enumerate(cfg.spectral)
        c = coefficients === nothing ? nothing : coefficients[s]
        for q in cfg.quantities, (w, λ) in enumerate(set.wavelengths)
            prefix, word = OPTICAL_QUANTITIES[q]
            k = if c === nothing
                zeros(nb)
            elseif q == :sct
                c.k_scat[w, :]
            elseif q == :abs
                c.k_abs[w, :]
            else
                c.k_scat[w, :] .+ c.k_abs[w, :]
            end
            lab = wavelength_label(λ)
            nm = string(round(Int, λ * 1e9), " nm")
            for (is3d, dimtag, depth) in ((true, "3d", "AOD"), (false, "2d", "AOT"))
                (is3d ? cfg.layers : cfg.columns) || continue
                stem = string(prefix, dimtag, "_", lab, "_")
                long = string(word, " ", depth, " ", nm)
                if cfg.bin_variables
                    for (b, bin) in enumerate(bins)
                        kb = zeros(nb)
                        kb[b] = k[b]
                        push!(
                            specs,
                            OpticsVariableSpec(
                                stem * bin.name, long * " " * bin.name, "1", :tau, is3d, kb
                            ),
                        )
                    end
                end
                if cfg.class_variables
                    for (class, n) in classes
                        kc = [b <= n ? k[b] : 0.0 for b in 1:nb]
                        push!(
                            specs,
                            OpticsVariableSpec(
                                stem * class, class * " " * long, "1", :tau, is3d, kc
                            ),
                        )
                    end
                end
            end
        end
    end
    return specs
end

"""
    validate_optics_config(cfg::DustOpticsConfig) -> DustOpticsConfig

Checks the configuration and throws an `ArgumentError` describing the first problem
found: `cfg.stats` with [`validate_config`](@ref) (without input groups, with the
output variable names and the surface pressure of the pressure levels); at least one
bin; unique bin names; unique spectral-set names; known quantities; non-empty input
directory and name; unique output variable names. Returns `cfg`.

# Example

```julia
validate_optics_config(cfg);               # throws if cfg is inconsistent
```
"""
function validate_optics_config(cfg::DustOpticsConfig)
    isempty(cfg.bins) && throw(ArgumentError("no dust bins given"))
    allunique(b.name for b in cfg.bins) || throw(ArgumentError("bin names are not unique"))
    allunique(s.name for s in cfg.spectral) ||
        throw(ArgumentError("names of the spectral sets are not unique"))
    for q in cfg.quantities
        haskey(OPTICAL_QUANTITIES, q) || throw(
            ArgumentError(
                "unknown optical quantity :$q (known: $(collect(keys(OPTICAL_QUANTITIES))))",
            ),
        )
    end
    allunique(cfg.quantities) || throw(ArgumentError("quantities are not unique"))
    isempty(cfg.input_dir) && throw(ArgumentError("input_dir must not be empty"))
    isempty(cfg.input_name) && throw(ArgumentError("input_name must not be empty"))
    isempty(cfg.stats.groups) ||
        @warn "DustOpticsConfig: stats.groups are not used by the optical post-processing"
    isempty(cfg.stats.derived) ||
        @warn "DustOpticsConfig: stats.derived are not used by the optical post-processing"
    names = optics_variable_names(cfg)
    allunique(names) || throw(ArgumentError("output variable names are not unique"))
    pl = cfg.stats.pressure_levels
    checked = pl === nothing ? names : [names; pl.surface_pressure]
    validate_config(cfg.stats; names=checked, require_groups=false)
    return cfg
end
