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
