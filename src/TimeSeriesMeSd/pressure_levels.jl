# Interpolation from the hybrid σ–pressure model layers to pressure levels.

"""
    layer_pressure(ps, plm_l, l, vc::VerticalCoordinate) -> Float64

Mid-layer pressure of model layer `l` in a column with surface pressure `ps`.

For the σ layers (`l < vc.ls1`) the pressure follows the surface pressure,
`p = p_top + σ (ps − p_top)` with `σ = (plm_l − p_top) / (p_sf − p_top)`, where `plm_l` is
the reference mid-layer pressure of the layer (`plm` in the input files) and `p_sf`,
`p_top` are `vc.psf` and `vc.ptop`. For the pure pressure layers (`l ≥ vc.ls1`) the
pressure is `plm_l`. All pressures are in the same unit (hPa).

# Arguments

  - `ps::Real`: surface pressure of the column.
  - `plm_l::Real`: reference mid-layer pressure of layer `l`.
  - `l::Integer`: layer index (1 = lowest layer).
  - `vc::VerticalCoordinate`: parameters of the vertical coordinate.

# Example

```julia
julia> vc = VerticalCoordinate(31, 6.86, 0.2);

julia> layer_pressure(6.86, 5.0, 1, vc)      # reference column: plm itself
5.0

julia> layer_pressure(8.0, 5.0, 1, vc) ≈ 0.2 + 4.8 * 7.8 / 6.66
true
```
"""
@inline function layer_pressure(ps::Real, plm_l::Real, l::Integer, vc::VerticalCoordinate)
    l >= vc.ls1 && return Float64(plm_l)
    sigma = (plm_l - vc.ptop) / (vc.psf - vc.ptop)
    return vc.ptop + sigma * (ps - vc.ptop)
end

"""
    interpolate_to_pressure(x, ps, plm, vc, plev) -> Array{Float32,4}

Interpolates the monthly 3D field `x` from the model layers to the pressure levels `plev`.

For every column and month, the mid-layer pressures `p₁ > p₂ > … > p_L` of the model
layers are computed from the surface pressure `ps` of that column and month with
[`layer_pressure`](@ref). For a pressure level `P`:

  - `P > ps` (below the local surface): missing (`NaN32`);
  - `ps ≥ P ≥ p₁` (between the surface and the mid-point of the lowest layer): value of
    the lowest layer;
  - `p_l ≥ P ≥ p_{l+1}`: linear interpolation in `ln p` between layers `l` and `l + 1`,
    `x = x_l + (ln P − ln p_l) / (ln p_{l+1} − ln p_l) · (x_{l+1} − x_l)`;
  - `P < p_L` (above the mid-point of the uppermost layer): missing.

The result is missing where a value of the two bracketing layers or the surface pressure
is missing. The loop over the columns is multi-threaded.

# Arguments

  - `x::AbstractArray{<:Real,4}`: field `(lon, lat, layer, time)`, layer 1 at the bottom;
    missing values `NaN`.
  - `ps::AbstractArray{<:Real,3}`: surface pressure `(lon, lat, time)`, same unit as
    `plm` and `plev` (hPa).
  - `plm::AbstractVector{<:Real}`: reference mid-layer pressures of the layers.
  - `vc::VerticalCoordinate`: parameters of the vertical coordinate.
  - `plev::AbstractVector{<:Real}`: target pressure levels.

# Returns

`Array{Float32,4}` of size `(lon, lat, length(plev), time)`.

# Example

```julia
julia> vc = VerticalCoordinate(3, 10.0, 1.0);   # layers 1, 2 σ; layer 3 pressure

julia> plm = [9.0, 5.0, 0.5];

julia> x = reshape(Float32[log(9.0), log(5.0), log(0.5)], 1, 1, 3, 1);  # x = ln p

julia> ps = fill(10.0f0, 1, 1, 1);

julia> y = interpolate_to_pressure(x, ps, plm, vc, [12.0, 9.5, 7.0, 2.0, 0.1]);

julia> isnan(y[1]), y[2] ≈ log(9.0), y[3] ≈ log(7.0), y[4] ≈ log(2.0), isnan(y[5])
(true, true, true, true, true)
```
"""
function interpolate_to_pressure(
    x::AbstractArray{<:Real,4},
    ps::AbstractArray{<:Real,3},
    plm::AbstractVector{<:Real},
    vc::VerticalCoordinate,
    plev::AbstractVector{<:Real},
)
    nlon, nlat, nlev, nt = size(x)
    size(ps) == (nlon, nlat, nt) || throw(
        DimensionMismatch(
            "surface pressure has size $(size(ps)), expected $((nlon, nlat, nt))"
        ),
    )
    length(plm) == nlev ||
        throw(DimensionMismatch("plm has $(length(plm)) levels, field has $nlev"))
    all(>(0), plm) || throw(ArgumentError("plm must be positive"))
    issorted(plm; rev=true) || throw(
        ArgumentError("plm must decrease with the layer index (layer 1 at the bottom)")
    )
    np = length(plev)
    target = Float64.(plev)
    logtarget = log.(target)
    out = Array{Float32}(undef, nlon, nlat, np, nt)
    @threads for jt in 1:(nlat * nt)
        j = mod1(jt, nlat)
        t = cld(jt, nlat)
        logp = Vector{Float64}(undef, nlev)
        for i in 1:nlon
            psc = Float64(ps[i, j, t])
            if isnan(psc)
                @inbounds for k in 1:np
                    out[i, j, k, t] = NaN32
                end
                continue
            end
            @inbounds for l in 1:nlev
                logp[l] = log(layer_pressure(psc, plm[l], l, vc))
            end
            @inbounds for k in 1:np
                P = target[k]
                lP = logtarget[k]
                value = if P > psc || lP < logp[nlev]
                    NaN
                elseif lP >= logp[1]
                    Float64(x[i, j, 1, t])
                else
                    l = 1
                    while logp[l + 1] > lP
                        l += 1
                    end
                    w = (lP - logp[l]) / (logp[l + 1] - logp[l])
                    x1 = Float64(x[i, j, l, t])
                    x2 = Float64(x[i, j, l + 1, t])
                    x1 + w * (x2 - x1)
                end
                out[i, j, k, t] = Float32(value)
            end
        end
    end
    return out
end

"""
    pressure_level_bounds(levels) -> Matrix{Float64}

Pressures at the edges of the layers represented by the pressure levels `levels`, as
CF cell bounds of size `(2, length(levels))`: `bounds[1, k]` is the lower edge (higher
pressure) and `bounds[2, k]` the upper edge (lower pressure) of the layer of level
`levels[k]`.

With the levels sorted by decreasing pressure, `p₁ > p₂ > … > p_n`, the edge between two
neighboring levels is their geometric mean, i.e. the midpoint in `ln p` (the coordinate
in which [`interpolate_to_pressure`](@ref) interpolates linearly),

    p_{k+1/2} = √(p_k p_{k+1}),

and the outer edges lie half a level spacing in `ln p` beyond the outermost levels,
`p_{1/2} = p₁ √(p₁/p₂)` and `p_{n+1/2} = p_n √(p_n/p_{n−1})`. Adjacent layers share their
edges, all edges are positive, and every level lies inside its layer. The levels may be
given in any order; the bounds refer to `levels[k]` in the given order. At least two
levels are needed.

# Arguments

  - `levels::AbstractVector{<:Real}`: positive, unique pressure levels.

# Returns

`Matrix{Float64}` of size `(2, length(levels))`.

# Example

```julia
julia> b = pressure_level_bounds([4.0, 1.0, 0.25]);

julia> b[:, 2]                     # layer of 1 hPa: edges √(4·1) and √(1·0.25)
2-element Vector{Float64}:
 2.0
 0.5

julia> b[1, 1], b[2, 3]            # outer edges: 4·√4 and 0.25·√0.25
(8.0, 0.125)
```
"""
function pressure_level_bounds(levels::AbstractVector{<:Real})
    check_pressure_levels(levels)
    n = length(levels)
    n >= 2 || throw(ArgumentError("pressure-level bounds need at least two levels"))
    order = sortperm(levels; rev=true)
    p = Float64.(levels[order])                # decreasing pressure
    edges = Vector{Float64}(undef, n + 1)      # edges[k] below, edges[k + 1] above p[k]
    for k in 1:(n - 1)
        edges[k + 1] = sqrt(p[k] * p[k + 1])
    end
    edges[1] = p[1] * sqrt(p[1] / p[2])
    edges[n + 1] = p[n] * sqrt(p[n] / p[n - 1])
    bounds = Matrix{Float64}(undef, 2, n)
    for (k, i) in enumerate(order)
        bounds[1, i] = edges[k]
        bounds[2, i] = edges[k + 1]
    end
    return bounds
end

"""
    model_level_bounds(plm, ple, psf) -> Matrix{Float64}

Reference pressures at the edges of the model layers as CF cell bounds of size
`(2, L)`: `bounds[1, l]` is the lower edge of layer `l` (`psf` for the lowest layer, the
upper edge `ple[l - 1]` of the layer below otherwise), `bounds[2, l] = ple[l]` its upper
edge. These are the edge pressures of the hybrid σ–pressure layers in a column with
surface pressure `psf` (`VerticalCoordinate.psf`); `ple` is the variable of the same name
in the `scaleacc` files, which also defines σ at the layer edges.

Throws an `ArgumentError` unless `plm` and `ple` have the same length and every
reference mid-layer pressure lies inside its layer, `lower edge > plm[l] > ple[l] ≥ 0`.

# Arguments

  - `plm::AbstractVector{<:Real}`: reference mid-layer pressures (layer 1 at the bottom).
  - `ple::AbstractVector{<:Real}`: reference upper-edge pressures of the layers.
  - `psf::Real`: reference surface pressure.

# Example

```julia
julia> model_level_bounds([5.0, 1.0], [3.0, 0.2], 6.86)
2×2 Matrix{Float64}:
 6.86  3.0
 3.0   0.2
```
"""
function model_level_bounds(
    plm::AbstractVector{<:Real}, ple::AbstractVector{<:Real}, psf::Real
)
    n = length(plm)
    length(ple) == n || throw(ArgumentError("plm has $n levels, ple has $(length(ple))"))
    bounds = Matrix{Float64}(undef, 2, n)
    for l in 1:n
        lower = l == 1 ? Float64(psf) : Float64(ple[l - 1])
        upper = Float64(ple[l])
        lower > plm[l] > upper >= 0 || throw(
            ArgumentError(
                "layer $l: plm = $(plm[l]) is not between the edges $lower and $upper"
            ),
        )
        bounds[1, l] = lower
        bounds[2, l] = upper
    end
    return bounds
end

"""
    output_pressure_bounds(cfg::TimeSeriesConfig, plm=nothing, ple=nothing)
        -> Union{Nothing,Matrix{Float64}}

Pressures (hPa) at the edges of the layers of the output pressure levels
[`output_pressure_levels`](@ref)`(cfg, plm)`, as CF cell bounds of size
`(2, number of levels)`, written as `plev_bnds`:

  - for `levels = :model` and the upper-edge pressures `ple` of the input files given
    (e.g. `rd.grid.ple`): the reference edge pressures of the model layers,
    [`model_level_bounds`](@ref)`(plm, ple, cfg.vertical.psf)`;
  - otherwise: the midpoints in `ln p` between the levels,
    [`pressure_level_bounds`](@ref).

Returns `nothing` for a single level without `ple` (no edges can be derived).

# Arguments

  - `cfg::TimeSeriesConfig`: configuration with `pressure_levels`.
  - `plm = nothing`, `ple = nothing`: reference mid-layer and upper-edge pressures of the
    input files (hPa).

# Example

```julia
julia> cfgm = TimeSeriesConfig(; cfg_kwargs..., vertical=VerticalCoordinate(31, 6.86, 0.2),
                               pressure_levels=PressureLevels(:model));

julia> output_pressure_bounds(cfgm, Float32[5, 1], Float32[3, 0.2])[:, 1]
2-element Vector{Float64}:
 6.86
 3.0
```
"""
function output_pressure_bounds(cfg::TimeSeriesConfig, plm=nothing, ple=nothing)
    levels = output_pressure_levels(cfg, plm)
    if cfg.pressure_levels.levels === :model && ple !== nothing
        cfg.vertical === nothing &&
            throw(ArgumentError("pressure_levels requires the vertical coordinate"))
        return model_level_bounds(levels, ple, cfg.vertical.psf)
    end
    length(levels) >= 2 || return nothing
    return pressure_level_bounds(levels)
end

"""
    pressure_flags(cfg::TimeSeriesConfig) -> Dict{String,Bool}

`to_pressure` flag of every configured variable (input and derived), by name.

# Example

```julia
julia> pressure_flags(cfg)["airmass"]
false
```
"""
function pressure_flags(cfg::TimeSeriesConfig)
    flags = Dict{String,Bool}(
        v.name => v.to_pressure for g in cfg.groups for v in g.variables
    )
    for d in cfg.derived
        flags[d.name] = d.to_pressure
    end
    return flags
end

"""
    PRESSURE_INTERPOLATION_NOTE

Text of the global attribute `Vertical_interpolation` of output files with variables on
pressure levels.

# Example

```julia
julia> occursin("ln(p)", PRESSURE_INTERPOLATION_NOTE)
true
```
"""
const PRESSURE_INTERPOLATION_NOTE =
    "Monthly means interpolated from the hybrid sigma-pressure model layers to " *
    "pressure levels, linear in ln(p) between layer mid-points; levels below the " *
    "surface and above the uppermost layer mid-point are missing; between the surface " *
    "and the lowest layer mid-point the value of the lowest layer is used."
