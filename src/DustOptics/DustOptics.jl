## Author: Jan P. Perlwitz, Climate, Aerosol, and Pollution Research, LLC
##         (jan.p.perlwitz@caprllc.com)
## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-11-04 to 2026-11-08

"""
    DustOptics

Submodule of `PlanetPostStats` for the optical post-processing of the dust aerosol of
ModelE/ROCKE-3D simulations. Its input are the monthly time series of the dust tracer
mass mixing ratios and of the air mass of the model layers written by
[`TimeSeriesMeSd`](@ref PlanetPostStats.TimeSeriesMeSd) (e.g. `*.ijl.mts.nc` of the
dust tracers `Clay`, `Silt1`, …). For every run and month it computes

  - the dust mass of every size bin of the radiative calculations in every layer,
    `m_b = f_b q_s m_air`, from the mass mixing ratio `q_s` of its source tracer, the
    fraction `f_b` of the tracer in the bin (e.g. the split of the clay tracer into four
    sub-bins), and the air mass `m_air` of the layer ([`DustBin`](@ref),
    [`DustBinMasses`](@ref));
  - the effective radius of the size distribution in the layers, of the column dust
    load, and its mass-weighted column average ([`effective_radius`](@ref),
    [`column_mean_effective_radius`](@ref));
  - the scattering, absorption, and extinction optical depths `τ = k m_b` of every bin
    and of the cumulative size classes (`SumClay1Dust`, …, `PM2Dust`, …, `PM32Dust`) in
    the layers and of the column, at selected diagnostic wavelengths, with the mass
    coefficients `k` interpolated from the optical lookup tables of the radiation code
    to the wavelengths and the effective radii of the bins ([`SpectralTable`](@ref),
    [`SpectralSet`](@ref), [`interpolate_coefficients`](@ref));
  - the density-scaled opacity of every bin and size class in the layers, the optical
    depth of a layer divided by its air mass, `τ / m_air` (m² kg⁻¹), which, unlike the
    layer optical depth, can be interpolated to pressure levels
    ([`density_scaled_opacity`](@ref)).

All these variables are then written as time series, multi-year means, and standard
deviations of gridded, zonal-mean, and global-mean fields with the pipeline of
`TimeSeriesMeSd` (`process_runs`): one run at a time, and one variable at a time, since
every variable is a `LazyVariable` computed from the bin masses only when it is written.
The memory needed is that of the source tracer masses of one run plus a few fields.

All settings are collected in a [`DustOpticsConfig`](@ref); the entry point is
[`run_dust_optics`](@ref).

# Example

```julia
using PlanetPostStats.TimeSeriesMeSd
using PlanetPostStats.DustOptics

stats = TimeSeriesConfig(; simpath="/path/to/simulations", runs=["MarsDustRun1"],
    sim_period=(YearMonth(1, 1), YearMonth(20, 12)),
    stat_period=(YearMonth(3, 1), YearMonth(20, 12)), calendar=MARS_CALENDAR,
    seasons=MARS_SEASONS, vertical=VerticalCoordinate(31, 6.86, 0.2),
    output_dir="Traceroptics3d", output_name="traceroptics3d")
bins = [DustBin("Clay1"; source="Clay", radius=0.132e-6, density=2.5e3),
        DustBin("Clay2"; source="Clay", radius=0.23e-6, density=2.5e3),
        DustBin("Clay3"; source="Clay", radius=0.416e-6, density=2.5e3),
        DustBin("Clay4"; source="Clay", radius=0.766e-6, density=2.5e3),
        DustBin("Silt1"; radius=1.386e-6, density=2.65e3)]
sw = SpectralSet("sw", read_spectral_table("aer_sw_dust.nc", "sw_bands.csv"),
                 [0.55e-6, 0.89e-6])
cfg = DustOpticsConfig(; stats=stats, bins=bins, spectral=[sw],
    fractions=bin_fractions(bins, "Clay" => [0.009, 0.081, 0.234, 0.676]),
    input_dir="Tracermass3d", input_name="tracermass3d")
results = run_dust_optics(cfg)
annual_summary(results)
```
"""
module DustOptics

using CSV: CSV
using NCDatasets: NCDatasets, NCDataset, dimnames
using OrderedCollections: OrderedDict
using PlanetTools.DerivePhysicsVariables: DerivePhysicsVariables
using PlanetTools.PlanetAerosolUtils: interpolate_extinction_coeff_log
using PlanetTools.ReadWritePlanet: writenetcdf
using PlanetTools.StatisticUtils: calc_weighted_mean

using ..TimeSeriesMeSd: TimeSeriesMeSd
using ..TimeSeriesMeSd:
    AnyVariable,
    GridInfo,
    LazyVariable,
    RunData,
    TimeSeriesConfig,
    VariableData,
    YearMonth,
    process_runs,
    validate_config,
    yearmonth_range

include("bins.jl")
include("spectral.jl")
include("config.jl")
include("input.jl")
include("driver.jl")

# Size bins and array kernels
export DustBin, bin_sources, bin_fractions, optics_classes
export DustBinMasses, combine_bins, column_masses
export effective_radius, column_mean_effective_radius, density_scaled_opacity

# Optical lookup tables
export SpectralTable, read_band_bounds, read_spectral_table, wavelength_label
export SpectralSet, OpticalCoefficients, interpolate_coefficients, write_coefficients

# Configuration, input, and driver
export OPTICAL_QUANTITIES, DustOpticsConfig, run_fractions, validate_optics_config
export optics_variable_names, input_file, load_optics_run, optics_variables
export run_dust_optics, optical_coefficients

end # module DustOptics
