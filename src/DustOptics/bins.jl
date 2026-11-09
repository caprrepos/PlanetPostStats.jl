# Dust size bins of the radiative calculations, their masses, and the array kernels for
# effective radii and optical depths.

"""
    DustBin(name; source=name, radius, density)

One size bin of the dust aerosol as used in the optical calculations.

The mass of the bin is a fixed fraction of the mass of one transported dust tracer, its
`source`: several bins may share a source, e.g. the four sub-bins `Clay1`, …, `Clay4`
into which the clay tracer `Clay` of ModelE/ROCKE-3D is split for the radiative transfer,
while every silt bin is its own source (`Silt1` from `Silt1`, …). The fractions are given
per run ([`DustOpticsConfig`](@ref)`.fractions`, [`bin_fractions`](@ref)).

# Arguments

  - `name::AbstractString`: name of the bin, used in the output variable names
    (e.g. `ext3d_550nm_Clay1`).

# Keyword arguments

  - `source::AbstractString = name`: name of the dust tracer (mass mixing ratio) in the
    input files that contains the mass of this bin.
  - `radius::Real`: effective radius of the bin (m), used for the effective radius of the
    size distribution and to interpolate the optical properties to the bin.
  - `density::Real`: particle density (kg m⁻³).

# Example

```julia
julia> b = DustBin("Clay1"; source="Clay", radius=0.132e-6, density=2.5e3);

julia> b.source, b.radius
("Clay", 1.32e-7)
```
"""
struct DustBin
    name::String
    source::String
    radius::Float64
    density::Float64
    function DustBin(
        name::AbstractString; source::AbstractString=name, radius::Real, density::Real
    )
        isempty(name) && throw(ArgumentError("bin name must not be empty"))
        isempty(source) && throw(ArgumentError("bin $name: source must not be empty"))
        (isfinite(radius) && radius > 0) ||
            throw(ArgumentError("bin $name: radius must be positive"))
        (isfinite(density) && density > 0) ||
            throw(ArgumentError("bin $name: density must be positive"))
        return new(String(name), String(source), Float64(radius), Float64(density))
    end
end

"""
    bin_sources(bins) -> Vector{String}

Names of the dust tracers that are the sources of `bins`, unique, in order of first
appearance.

# Example

```julia
julia> bins = [DustBin("Clay1"; source="Clay", radius=1e-7, density=2500),
               DustBin("Clay2"; source="Clay", radius=2e-7, density=2500),
               DustBin("Silt1"; radius=1.4e-6, density=2650)];

julia> bin_sources(bins)
2-element Vector{String}:
 "Clay"
 "Silt1"
```
"""
bin_sources(bins::AbstractVector{DustBin}) = unique(b.source for b in bins)

"""
    bin_fractions(bins, splits::Pair...) -> Vector{Float64}

Fractions of the source tracer mass in each of the `bins`, for
[`DustOpticsConfig`](@ref)`.fractions`. Every pair `source => fractions` gives the
fractions of the bins with that source, in the order in which these bins appear in
`bins`; all other bins get the fraction 1 (the bin is the whole tracer).

# Arguments

  - `bins::AbstractVector{DustBin}`: the size bins.
  - `splits`: pairs `source::String => fractions::AbstractVector{<:Real}`.

# Example

```julia
julia> bins = [DustBin("Clay1"; source="Clay", radius=1e-7, density=2500),
               DustBin("Clay2"; source="Clay", radius=2e-7, density=2500),
               DustBin("Silt1"; radius=1.4e-6, density=2650)];

julia> bin_fractions(bins, "Clay" => [0.25, 0.75])
3-element Vector{Float64}:
 0.25
 0.75
 1.0
```
"""
function bin_fractions(bins::AbstractVector{DustBin}, splits::Pair...)
    f = ones(length(bins))
    for (source, fractions) in splits
        idx = findall(b -> b.source == source, bins)
        isempty(idx) && throw(ArgumentError("no bin has the source $source"))
        length(fractions) == length(idx) || throw(
            ArgumentError(
                "$(length(fractions)) fractions given for the $(length(idx)) bins of " *
                "source $source",
            ),
        )
        f[idx] .= fractions
    end
    return f
end

"""
    optics_classes(bins) -> Vector{Pair{String,Int}}

Cumulative size classes of the optical variables: pairs `class => n`, where the class
`class` is the sum over the first `n` bins of `bins`.

The classes and their names are those of
`PlanetTools.DerivePhysicsVariables.accumulate_dust_optics_modele`, for the bins
`Clay1`, `Clay2`, `Clay3`, `Clay4`, `Silt1`, …, `Silt5`:

| class          | sum of the bins           |
|:---------------|:--------------------------|
| `SumClay1Dust` | `Clay1`                   |
| `SumClay2Dust` | `Clay1 + Clay2`           |
| `SumClay3Dust` | `Clay1 + Clay2 + Clay3`   |
| `PM2Dust`      | `Clay1 + … + Clay4`       |
| `PM4Dust`      | … `+ Silt1`               |
| `PM8Dust`      | … `+ Silt2`               |
| `PM16Dust`     | … `+ Silt3`               |
| `PM32Dust`     | … `+ Silt4`               |
| `PM64Dust`     | … `+ Silt5`               |

As in PlanetTools, the accumulation stops at the first bin of `bins` that does not match
this sequence, so that only complete classes are created; for bins with other names the
result is empty.

# Example

```julia
julia> names = ["Clay1", "Clay2", "Clay3", "Clay4", "Silt1"];

julia> bins = [DustBin(n; radius=1e-6, density=2500) for n in names];

julia> optics_classes(bins)
5-element Vector{Pair{String, Int64}}:
 "SumClay1Dust" => 1
 "SumClay2Dust" => 2
 "SumClay3Dust" => 3
      "PM2Dust" => 4
      "PM4Dust" => 5
```
"""
function optics_classes(bins::AbstractVector{DustBin})
    seq = DerivePhysicsVariables.OPTICS_BINS
    classes = DerivePhysicsVariables.OPTICS_CLASSES
    out = Pair{String,Int}[]
    for (n, b) in enumerate(bins)
        (n <= length(seq) && b.name == seq[n]) || break
        push!(out, String(classes[n]) => n)
    end
    return out
end

"""
    DustBinMasses(sources, masses, source_index, fractions)

Dust masses per unit area in the model layers of one run, held once per source tracer:
the mass of bin `b` is `fractions[b] .* masses[source_index[b]]`.

# Fields

  - `sources::Vector{String}`: names of the source tracers.
  - `masses::Vector{Array{Float32,4}}`: mass per unit area of each source tracer
    (kg m⁻²), `(lon, lat, level, time)`, i.e. the mass mixing ratio times the air mass of
    the layer; missing values `NaN32`.
  - `source_index::Vector{Int}`: index into `sources` of the source of every bin.
  - `fractions::Vector{Float64}`: fraction of the source mass in every bin.

Linear combinations of the bin masses, `Σ_b c_b m_b`, are computed with
[`combine_bins`](@ref), so that the per-bin fields never have to be stored.

# Example

```julia
m = DustBinMasses(["Clay"], [rand(Float32, 4, 3, 2, 12)], [1, 1], [0.4, 0.6])
size(combine_bins(m, [1.0, 1.0]))      # (4, 3, 2, 12): total mass of the two bins
```
"""
struct DustBinMasses
    sources::Vector{String}
    masses::Vector{Array{Float32,4}}
    source_index::Vector{Int}
    fractions::Vector{Float64}
    function DustBinMasses(sources, masses, source_index, fractions)
        length(sources) == length(masses) ||
            throw(ArgumentError("one mass array per source needed"))
        length(source_index) == length(fractions) ||
            throw(ArgumentError("source_index and fractions differ in length"))
        all(i -> 1 <= i <= length(sources), source_index) ||
            throw(ArgumentError("source index out of range"))
        isempty(masses) ||
            all(m -> size(m) == size(masses[1]), masses) ||
            throw(DimensionMismatch("mass arrays of the sources differ in size"))
        return new(
            collect(String, sources),
            collect(Array{Float32,4}, masses),
            collect(Int, source_index),
            collect(Float64, fractions),
        )
    end
end

"""
    combine_bins(m::DustBinMasses, coefficients) -> Array{Float32,4}

The linear combination `Σ_b coefficients[b] · fractions[b] · masses[source of b]` of the
bin masses, computed per source tracer: the coefficients of all bins of one source are
first summed, so every source array is read only once. Bins with a zero coefficient do
not contribute, i.e. a missing value of their source does not propagate. Accumulation in
`Float32`.

This one operation gives all fields of the optical post-processing:

  - optical depth of bin `b` or of a class of bins: coefficients `k_b` (mass extinction,
    scattering, or absorption coefficient, m² kg⁻¹) for the bins of the class, `0`
    otherwise;
  - total dust mass of the layer: coefficients `1`;
  - numerator and denominator of the effective radius: `1/ρ_b` and `1/(ρ_b r_b)`
    (see [`effective_radius`](@ref)).

# Arguments

  - `m::DustBinMasses`: the bin masses.
  - `coefficients::AbstractVector{<:Real}`: one coefficient per bin.

# Example

```julia
m = DustBinMasses(["Clay", "Silt1"], [fill(1.0f0, 1, 1, 1, 1), fill(2.0f0, 1, 1, 1, 1)],
                  [1, 1, 2], [0.4, 0.6, 1.0])
combine_bins(m, [10.0, 20.0, 1.0])[1]      # 10 · 0.4 · 1 + 20 · 0.6 · 1 + 1 · 2 = 18
```
"""
function combine_bins(m::DustBinMasses, coefficients::AbstractVector{<:Real})
    length(coefficients) == length(m.fractions) ||
        throw(ArgumentError("one coefficient per bin needed"))
    isempty(m.masses) && throw(ArgumentError("no source masses"))
    c = zeros(length(m.sources))
    for (b, s) in enumerate(m.source_index)
        c[s] += coefficients[b] * m.fractions[b]
    end
    out = zeros(Float32, size(m.masses[1]))
    for (s, cs) in enumerate(c)
        iszero(cs) && continue
        out .+= Float32(cs) .* m.masses[s]
    end
    return out
end

"""
    column_masses(m::DustBinMasses) -> DustBinMasses

Column masses of the sources, i.e. the masses summed over the layers (dimension 3), as a
[`DustBinMasses`](@ref) with a vertical dimension of length 1. A missing value in any
layer makes the column missing. [`combine_bins`](@ref) of the result gives column
quantities, e.g. the column optical depth.

# Example

```julia
m = DustBinMasses(["Clay"], [ones(Float32, 4, 3, 2, 12)], [1], [1.0])
size(column_masses(m).masses[1])           # (4, 3, 1, 12)
```
"""
function column_masses(m::DustBinMasses)
    cols = [sum(x; dims=3) for x in m.masses]
    return DustBinMasses(m.sources, cols, m.source_index, m.fractions)
end

"""
    density_scaled_opacity(m::DustBinMasses, coefficients, airmass) -> Array{Float32,4}

Density-scaled opacity (m² kg⁻¹) of the layers: the optical depth of a layer divided by
the air mass of the layer per unit area,

    dso = τ / m_air = Σ_b k_b m_b / m_air = Σ_b k_b f_b q_s,

with the bin masses `m_b = f_b q_s m_air` of `m`, the mass coefficients
`coefficients` (`k_b`, m² kg⁻¹; zero for bins not included, see
[`combine_bins`](@ref)), and the air mass `airmass` (kg m⁻²) used to form `m`. It equals
the extinction (or scattering, absorption) coefficient per unit volume divided by the
air density, `β / ρ_air`, and the optical depth per unit pressure divided by the
gravitational acceleration, `(dτ/dp) / g`. Unlike the optical depth of a layer, it does
not depend on the thickness of the layer: it is proportional to the dust mass mixing
ratio and can be interpolated to pressure levels. The Mars Climate Sounder reports the
mass mixing ratio of dust from its retrieved opacity profiles in this form [Heavens et al.,
2011]; for the opacity per unit length itself see [`opacity_per_km`](@ref).

Cells with a missing value or a non-positive air mass are `NaN32`.

# Arguments

  - `m::DustBinMasses`: bin masses of the layers.
  - `coefficients::AbstractVector{<:Real}`: one coefficient per bin (m² kg⁻¹).
  - `airmass::AbstractArray{<:Real,4}`: air mass of the layers (kg m⁻²), same size as
    the arrays of `m`.

# References

  - Heavens, N. G., et al. (2011), The vertical distribution of dust in the Martian
    atmosphere during northern spring and summer: Observations by the Mars Climate
    Sounder and analysis of zonal average vertical dust profiles, *J. Geophys. Res.*,
    116, E04003, doi:10.1029/2010JE003691.

# Example

```julia
air = fill(2.0f0, 1, 1, 1, 1)
m = DustBinMasses(["Clay"], [0.5f0 .* air], [1], [1.0])   # q = 0.5 kg kg⁻¹
density_scaled_opacity(m, [3.0], air)[1]                   # 3 · 0.5 = 1.5 m² kg⁻¹
```
"""
function density_scaled_opacity(
    m::DustBinMasses, coefficients::AbstractVector{<:Real}, airmass::AbstractArray{<:Real,4}
)
    isempty(m.masses) ||
        size(airmass) == size(m.masses[1]) ||
        throw(
            DimensionMismatch(
                "airmass has size $(size(airmass)), masses $(size(m.masses[1]))"
            ),
        )
    tau = combine_bins(m, coefficients)
    @. tau = ifelse(airmass > 0, tau / airmass, NaN32)
    return tau
end

"""
    air_density(ps, temperature, plm, vc::VerticalCoordinate, gas_constant)
        -> Array{Float32,4}

Air density (kg m⁻³) at the mid-points of the model layers, from the ideal gas law,

    ρ = p / (R T),

with the mid-layer pressure `p` of the hybrid σ–pressure coordinate
(`TimeSeriesMeSd.layer_pressure`: `p = p_top + σ (p_s − p_top)` below layer `LS1`,
`p = plm` above), the temperature `T` of the layer, and the specific gas constant `R`
of the atmosphere.

# Arguments

  - `ps::AbstractArray{<:Real,3}`: surface pressure (hPa), `(lon, lat, time)`.
  - `temperature::AbstractArray{<:Real,4}`: layer temperature (K),
    `(lon, lat, level, time)`.
  - `plm::AbstractVector{<:Real}`: reference mid-layer pressures (hPa).
  - `vc::VerticalCoordinate`: parameters of the vertical coordinate.
  - `gas_constant::Real`: specific gas constant (J kg⁻¹ K⁻¹), e.g. 188.92 for CO₂.

Cells with a missing value or a non-positive temperature are `NaN32`.

# Example

```julia
vc = VerticalCoordinate(31, 6.86, 0.2)
ρ = air_density(fill(6.86f0, 1, 1, 1), fill(200.0f0, 1, 1, 1, 1), [5.0], vc, 188.92)
ρ[1]                                    # 500 Pa / (188.92 J kg⁻¹ K⁻¹ · 200 K) ≈ 0.0132
```
"""
function air_density(
    ps::AbstractArray{<:Real,3},
    temperature::AbstractArray{<:Real,4},
    plm::AbstractVector{<:Real},
    vc::VerticalCoordinate,
    gas_constant::Real,
)
    nlon, nlat, nlev, nt = size(temperature)
    size(ps) == (nlon, nlat, nt) || throw(
        DimensionMismatch(
            "surface pressure has size $(size(ps)), expected $((nlon, nlat, nt))"
        ),
    )
    length(plm) == nlev ||
        throw(DimensionMismatch("plm has $(length(plm)) levels, temperature $nlev"))
    gas_constant > 0 || throw(ArgumentError("gas_constant must be positive"))
    rho = Array{Float32}(undef, nlon, nlat, nlev, nt)
    @inbounds for t in 1:nt, l in 1:nlev, j in 1:nlat, i in 1:nlon
        T = Float64(temperature[i, j, l, t])
        p = 100 * layer_pressure(Float64(ps[i, j, t]), plm[l], l, vc)   # Pa
        rho[i, j, l, t] = T > 0 && isfinite(p) ? Float32(p / (gas_constant * T)) : NaN32
    end
    return rho
end

"""
    opacity_per_km(m::DustBinMasses, coefficients, airmass, density) -> Array{Float32,4}

Opacity per kilometer (km⁻¹) of the layers, the extinction (or scattering, absorption)
coefficient per unit volume,

    β = dso · ρ_air · 1000 m km⁻¹,     dso = Σ_b k_b m_b / m_air,

i.e. the [`density_scaled_opacity`](@ref) times the air density `density` (kg m⁻³, e.g.
[`air_density`](@ref)). This is the quantity of the dust opacity profiles retrieved by
the Mars Climate Sounder at pressure levels [Heavens et al., 2011]. It does not depend
on the layer thickness and can be interpolated to pressure levels.

# Arguments

  - `m::DustBinMasses`: bin masses of the layers.
  - `coefficients::AbstractVector{<:Real}`: one coefficient per bin (m² kg⁻¹).
  - `airmass::AbstractArray{<:Real,4}`: air mass of the layers (kg m⁻²).
  - `density::AbstractArray{<:Real,4}`: air density of the layers (kg m⁻³).

# References

  - Heavens, N. G., et al. (2011), The vertical distribution of dust in the Martian
    atmosphere during northern spring and summer: Observations by the Mars Climate
    Sounder and analysis of zonal average vertical dust profiles, *J. Geophys. Res.*,
    116, E04003, doi:10.1029/2010JE003691.

# Example

```julia
air = fill(2.0f0, 1, 1, 1, 1)
m = DustBinMasses(["Clay"], [0.5f0 .* air], [1], [1.0])
opacity_per_km(m, [3.0], air, fill(0.01f0, 1, 1, 1, 1))[1]   # 1.5 · 0.01 · 1000 = 15 km⁻¹
```
"""
function opacity_per_km(
    m::DustBinMasses,
    coefficients::AbstractVector{<:Real},
    airmass::AbstractArray{<:Real,4},
    density::AbstractArray{<:Real,4},
)
    size(density) == size(airmass) || throw(
        DimensionMismatch("density has size $(size(density)), airmass $(size(airmass))")
    )
    beta = density_scaled_opacity(m, coefficients, airmass)
    @. beta *= 1000.0f0 * density
    return beta
end

"""
    effective_radius(m::DustBinMasses, bins) -> Array{Float32,4}

Effective radius (m) of the size distribution given by the bin masses `m` of the size
bins `bins`, in every grid cell and month,

    r_eff = Σ_b V_b / Σ_b (V_b / r_b),     V_b = m_b / ρ_b,

the area-weighted mean radius of spherical particles [Hansen and Travis, 1974], with the
effective radius `r_b` and the particle density `ρ_b` of every bin. This is the formula
of `PlanetTools.PlanetAerosolUtils.calc_size_integrated_effective_radius`, evaluated for
all grid cells at once with [`combine_bins`](@ref). Cells without dust (zero total mass)
and cells with a missing mass are `NaN32`.

Applied to [`column_masses`](@ref)`(m)`, the result is the effective radius of the
column dust load.

# Arguments

  - `m::DustBinMasses`: bin masses (layer or column).
  - `bins::AbstractVector{DustBin}`: the bins, in the order of `m.fractions`.

# References

  - Hansen, J. E., and L. D. Travis (1974), Light scattering in planetary atmospheres,
    *Space Sci. Rev.*, 16(4), 527–610, doi:10.1007/BF00168069.

# Example

```julia
bins = [DustBin("a"; radius=1e-7, density=2500), DustBin("b"; radius=1e-6, density=2500)]
m = DustBinMasses(["a", "b"], [ones(Float32, 1, 1, 1, 1), ones(Float32, 1, 1, 1, 1)],
                  [1, 2], [1.0, 1.0])
effective_radius(m, bins)[1]               # 2 / (1/1e-7 + 1/1e-6) ≈ 1.82e-7
```
"""
function effective_radius(m::DustBinMasses, bins::AbstractVector{DustBin})
    length(bins) == length(m.fractions) ||
        throw(ArgumentError("one bin per bin fraction needed"))
    volume = combine_bins(m, [1 / b.density for b in bins])
    surface = combine_bins(m, [1 / (b.density * b.radius) for b in bins])
    return volume ./ surface
end

"""
    column_mean_effective_radius(m::DustBinMasses, bins) -> Array{Float32,3}

Column average of the effective radius of the layers, weighted with the dust mass of
the layers,

    r̄_eff = Σ_l M_l r_eff,l / Σ_l M_l,     M_l = Σ_b m_b,l,

computed with `PlanetTools.StatisticUtils.calc_weighted_mean` over the vertical
dimension: layers without dust or with missing values are ignored, columns without any
valid layer are `NaN32`. Result `(lon, lat, time)`.

In general this differs from the effective radius of the column dust load,
[`effective_radius`](@ref)`(`[`column_masses`](@ref)`(m), bins)`, which weights the
layers with their particle cross sections instead of their masses.

# Example

```julia
r = column_mean_effective_radius(m, bins)
size(r)                                    # (nlon, nlat, nmonths)
```
"""
function column_mean_effective_radius(m::DustBinMasses, bins::AbstractVector{DustBin})
    reff = effective_radius(m, bins)
    mtot = combine_bins(m, ones(length(bins)))
    return dropdims(calc_weighted_mean(reff, mtot; dims=3); dims=3)
end
