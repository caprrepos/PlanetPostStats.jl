# Optical lookup tables of the radiation code and their interpolation to the diagnostic
# wavelengths and the size bins.

"""
    SpectralTable(band_bounds, radius, k_scat, k_abs, g_asym)

Optical properties of one aerosol component in the spectral bands of a radiation code,
tabulated for a set of particle effective radii (e.g. the aerosol optical-property files
of SOCRATES used by ROCKE-3D).

# Fields

  - `band_bounds::Vector{Tuple{Float64,Float64}}`: `(lower, upper)` wavelength bounds of
    the bands (m).
  - `radius::Vector{Float64}`: effective radii of the table (m).
  - `k_scat::Matrix{Float64}`: mass scattering coefficient (m² kg⁻¹), `(band, radius)`.
  - `k_abs::Matrix{Float64}`: mass absorption coefficient (m² kg⁻¹), `(band, radius)`.
  - `g_asym::Matrix{Float64}`: asymmetry parameter, `(band, radius)`.

Row `i` of the matrices belongs to band `i` of `band_bounds`. The constructor checks the
sizes and that the bounds and radii are positive.

See also [`read_spectral_table`](@ref), [`interpolate_coefficients`](@ref).

# Example

```julia
julia> t = SpectralTable([(0.4e-6, 0.6e-6), (0.6e-6, 1.0e-6)], [1e-7, 1e-6],
                         [1.0 2.0; 3.0 4.0], [0.1 0.2; 0.3 0.4], [0.6 0.7; 0.6 0.7]);

julia> size(t.k_scat)
(2, 2)
```
"""
struct SpectralTable
    band_bounds::Vector{Tuple{Float64,Float64}}
    radius::Vector{Float64}
    k_scat::Matrix{Float64}
    k_abs::Matrix{Float64}
    g_asym::Matrix{Float64}
    function SpectralTable(band_bounds, radius, k_scat, k_abs, g_asym)
        nb, nr = length(band_bounds), length(radius)
        for (name, x) in (("k_scat", k_scat), ("k_abs", k_abs), ("g_asym", g_asym))
            size(x) == (nb, nr) ||
                throw(DimensionMismatch("$name has size $(size(x)), expected ($nb, $nr)"))
        end
        all(b -> 0 < b[1] <= b[2], band_bounds) ||
            throw(ArgumentError("band bounds must satisfy 0 < lower ≤ upper"))
        all(>(0), radius) || throw(ArgumentError("radii must be positive"))
        bounds = [(Float64(b[1]), Float64(b[2])) for b in band_bounds]
        return new(
            bounds,
            Float64.(collect(radius)),
            Float64.(k_scat),
            Float64.(k_abs),
            Float64.(g_asym),
        )
    end
end

"""
    read_band_bounds(path; header_row=3, band_column="Band") -> Vector{Tuple{Float64,Float64}}

Reads the wavelength bounds of the spectral bands from the CSV file `path`, with CSV.jl.

The column names are in line `header_row` of the file (the lines above are skipped, e.g.
comments), followed by one line per band. Besides the column `band_column` with the band
number, the file must have exactly two numeric columns, the lower and upper wavelength
bound of the band (m); the order of the two does not matter. The bands are returned in
the order of the file, which must be the order of the bands in the optical-property
file (see [`read_spectral_table`](@ref)).

# Arguments

  - `path::AbstractString`: CSV file.

# Keyword arguments

  - `header_row::Integer = 3`: line with the column names.
  - `band_column::AbstractString = "Band"`: name of the column with the band numbers;
    ignored if absent.

# Example

```julia
bounds = read_band_bounds("sp_sw_bands.csv")
bounds[1]                                  # (lower, upper) of band 1 in m
```
"""
function read_band_bounds(
    path::AbstractString; header_row::Integer=3, band_column::AbstractString="Band"
)
    isfile(path) || throw(ArgumentError("band file not found: $path"))
    f = CSV.File(path; header=header_row)
    names = [String(n) for n in propertynames(f) if String(n) != band_column]
    length(names) == 2 || throw(
        ArgumentError(
            "band file $path must have two columns besides $band_column, found $names"
        ),
    )
    lo = Float64.(getproperty(f, Symbol(names[1])))
    hi = Float64.(getproperty(f, Symbol(names[2])))
    return [(min(a, b), max(a, b)) for (a, b) in zip(lo, hi)]
end

"""
    read_spectral_table(spectral_file, bands_file; header_row=3, humidity_index=1,
                        component_index=1) -> SpectralTable

Reads the optical properties of one aerosol component from the NetCDF file
`spectral_file` and the wavelength bounds of its bands from the CSV file `bands_file`
([`read_band_bounds`](@ref)).

Variables read from `spectral_file` (dimensions in Julia order):

  - `radius_eff_dry(radius, component)`: dry effective radius (m);
  - `k_scat`, `k_abs`, `g_asym(band, humidity, radius, component)`: mass scattering and
    absorption coefficients (m² kg⁻¹) and asymmetry parameter.

The table of humidity `humidity_index` and component `component_index` is used. The
number of bands of the two files must agree; the bands are paired in file order.

# Arguments

  - `spectral_file::AbstractString`: NetCDF file of optical properties.
  - `bands_file::AbstractString`: CSV file of band bounds.

# Keyword arguments

  - `header_row::Integer = 3`: header line of `bands_file`.
  - `humidity_index::Integer = 1`, `component_index::Integer = 1`: table to use.

# Example

```julia
sw = read_spectral_table("aer_sw_dust.nc", "sp_sw_bands.csv")
length(sw.band_bounds), length(sw.radius)
```
"""
function read_spectral_table(
    spectral_file::AbstractString,
    bands_file::AbstractString;
    header_row::Integer=3,
    humidity_index::Integer=1,
    component_index::Integer=1,
)
    isfile(spectral_file) ||
        throw(ArgumentError("optical-property file not found: $spectral_file"))
    bounds = read_band_bounds(bands_file; header_row=header_row)
    radius, ks, ka, ga = NCDataset(spectral_file, "r") do ds
        for v in ("radius_eff_dry", "k_scat", "k_abs", "g_asym")
            haskey(ds, v) || throw(ArgumentError("variable $v not found in $spectral_file"))
        end
        r = Array(ds["radius_eff_dry"])
        ndims(r) == 2 || throw(
            ArgumentError("radius_eff_dry in $spectral_file must be (radius, component)"),
        )
        tables = map(("k_scat", "k_abs", "g_asym")) do v
            x = Array(ds[v])
            ndims(x) == 4 || throw(
                ArgumentError(
                    "$v in $spectral_file must be (band, humidity, radius, component)"
                ),
            )
            return Float64.(coalesce.(x[:, humidity_index, :, component_index], NaN))
        end
        return (Float64.(coalesce.(r[:, component_index], NaN)), tables...)
    end
    size(ks, 1) == length(bounds) || throw(
        DimensionMismatch(
            "$spectral_file has $(size(ks, 1)) bands, $bands_file has $(length(bounds))"
        ),
    )
    return SpectralTable(bounds, radius, ks, ka, ga)
end

"""
    wavelength_label(wavelength) -> String

Label of a wavelength (m) in the names of the output variables: the wavelength in
nanometers, rounded to an integer, followed by `nm`.

# Example

```julia
julia> wavelength_label(0.55e-6), wavelength_label(21.6e-6)
("550nm", "21600nm")
```
"""
wavelength_label(wavelength::Real) = string(round(Int, wavelength * 1e9), "nm")

"""
    SpectralSet(name, table, wavelengths; coefficients_file=nothing)

Diagnostic wavelengths of one spectral range (e.g. shortwave or longwave) and the
optical lookup table from which their optical properties are interpolated.

# Arguments

  - `name::AbstractString`: name of the set (e.g. `"sw"`), used in the variable names of
    the coefficient file.
  - `table::SpectralTable`: optical properties in the bands of the radiation code.
  - `wavelengths::AbstractVector{<:Real}`: diagnostic wavelengths (m), positive, with
    unique [`wavelength_label`](@ref)s.

# Keyword arguments

  - `coefficients_file = nothing`: if a path is given, the coefficients interpolated to
    the wavelengths and the size bins are written to this NetCDF file
    ([`write_coefficients`](@ref)).

# Example

```julia
sw = SpectralSet("sw", read_spectral_table("aer_sw.nc", "sw_bands.csv"), [0.55e-6, 0.89e-6])
```
"""
struct SpectralSet
    name::String
    table::SpectralTable
    wavelengths::Vector{Float64}
    coefficients_file::Union{Nothing,String}
    function SpectralSet(
        name::AbstractString,
        table::SpectralTable,
        wavelengths::AbstractVector{<:Real};
        coefficients_file::Union{Nothing,AbstractString}=nothing,
    )
        isempty(wavelengths) && throw(ArgumentError("set $name: no wavelengths"))
        all(w -> isfinite(w) && w > 0, wavelengths) ||
            throw(ArgumentError("set $name: wavelengths must be positive"))
        allunique(wavelength_label.(wavelengths)) ||
            throw(ArgumentError("set $name: wavelength labels are not unique"))
        file = coefficients_file === nothing ? nothing : String(coefficients_file)
        return new(String(name), table, Float64.(collect(wavelengths)), file)
    end
end

"""
    OpticalCoefficients

Optical properties of the size bins at the diagnostic wavelengths of one
[`SpectralSet`](@ref), as returned by [`interpolate_coefficients`](@ref).

# Fields

  - `name::String`: name of the spectral set.
  - `wavelengths::Vector{Float64}`: diagnostic wavelengths (m).
  - `k_scat`, `k_abs::Matrix{Float64}`: mass scattering and absorption coefficients
    (m² kg⁻¹), `(wavelength, bin)`.
  - `g_asym::Matrix{Float64}`: asymmetry parameter, `(wavelength, bin)`.

# Example

```julia
c = interpolate_coefficients(sw, bins)
c.k_scat[1, :] .+ c.k_abs[1, :]           # mass extinction coefficients at λ₁
```
"""
struct OpticalCoefficients
    name::String
    wavelengths::Vector{Float64}
    k_scat::Matrix{Float64}
    k_abs::Matrix{Float64}
    g_asym::Matrix{Float64}
end

"""
    interpolate_coefficients(set::SpectralSet, bins; extrapolation=:linear)
        -> OpticalCoefficients

Interpolates the mass scattering and absorption coefficients and the asymmetry
parameter of `set.table` to the wavelengths `set.wavelengths` and the effective radii of
the size `bins`, with `PlanetTools.PlanetAerosolUtils.interpolate_extinction_coeff_log`:
linear in `log(value)` versus `log(wavelength)` and `log(radius)`, each band placed at
the geometric mean of its bounds.

# Arguments

  - `set::SpectralSet`: wavelengths and lookup table.
  - `bins::AbstractVector{DustBin}`: size bins.

# Keyword arguments

  - `extrapolation::Symbol = :linear`: treatment of wavelengths and radii outside the
    table, `:linear` (power laws), `:flat`, or `:error`.

# Example

```julia
c = interpolate_coefficients(sw, bins)
size(c.k_scat)                             # (length(sw.wavelengths), length(bins))
```
"""
function interpolate_coefficients(
    set::SpectralSet, bins::AbstractVector{DustBin}; extrapolation::Symbol=:linear
)
    t = set.table
    radii = [b.radius for b in bins]
    itp(x) = interpolate_extinction_coeff_log(
        t.band_bounds, t.radius, x, set.wavelengths, radii; extrapolation=extrapolation
    )
    return OpticalCoefficients(
        set.name, copy(set.wavelengths), itp(t.k_scat), itp(t.k_abs), itp(t.g_asym)
    )
end

"""
    write_coefficients(path, c::OpticalCoefficients, bins) -> String

Writes the interpolated optical properties `c` of the size `bins` to the NetCDF file
`path` (created or overwritten; parent directories are created), with
`PlanetTools.ReadWritePlanet.writenetcdf`: variables `diag_k_scat_<name>`,
`diag_k_abs_<name>` (m² kg⁻¹), and `diag_g_asym_<name>` on the dimensions `wavelength`
(m) and `radius_eff_dry` (m, the effective radii of the bins), where `<name>` is
`c.name`. Returns `path`.

# Example

```julia
write_coefficients("aer_sw_interpol.nc", interpolate_coefficients(sw, bins), bins)
```
"""
function write_coefficients(
    path::AbstractString, c::OpticalCoefficients, bins::AbstractVector{DustBin}
)
    mkpath(dirname(abspath(path)))
    isfile(path) && rm(path)
    dimkw = (;
        datadim1=Float32.(c.wavelengths),
        dim1attr=OrderedDict{String,Any}("units" => "m", "long_name" => "Wavelength"),
        datadim2=Float32[b.radius for b in bins],
        dim2attr=OrderedDict{String,Any}(
            "units" => "m", "long_name" => "Dry effective radius"
        ),
        dimnames=["wavelength", "radius_eff_dry"],
    )
    vars = (
        ("diag_k_scat_", c.k_scat, "m2 kg-1", "Mass scattering coefficient"),
        ("diag_k_abs_", c.k_abs, "m2 kg-1", "Mass absorption coefficient"),
        ("diag_g_asym_", c.g_asym, "1", "Asymmetry parameter"),
    )
    for (i, (prefix, data, units, long)) in enumerate(vars)
        writenetcdf(
            path;
            mode=i == 1 ? "c" : "a",
            data=Float32.(data),
            names=prefix * c.name,
            dataattr=OrderedDict{String,Any}("units" => units, "long_name" => long),
            dimkw...,
        )
    end
    return String(path)
end
