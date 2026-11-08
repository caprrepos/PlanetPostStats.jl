# Entry point of the optical post-processing.

"""
    run_dust_optics(cfg::DustOpticsConfig; runs=cfg.stats.runs) -> Vector{RunResult}

Main entry point of the optical post-processing of dust. Validates `cfg`
([`validate_optics_config`](@ref)) and

 1. interpolates the optical properties of every spectral set to its diagnostic
    wavelengths and to the effective radii of the size bins
    ([`interpolate_coefficients`](@ref)) and writes them to the `coefficients_file` of
    the set, if one is given ([`write_coefficients`](@ref));
 2. processes the runs one after the other with
    `TimeSeriesMeSd.process_runs`: reads the dust tracers and the air mass of the run
    ([`load_optics_run`](@ref)) and computes and writes every output variable, one at a
    time, with the time series and statistics of `cfg.stats` (output in
    `output_directory(cfg.stats, run)`).

Returns the `RunResult`s of the runs (written files, summary tables of the global means
of the 2D variables; see `annual_summary`, `seasonal_summary`).

# Arguments

  - `cfg::DustOpticsConfig`: configuration.

# Keyword arguments

  - `runs = cfg.stats.runs`: runs to process.

# Example

```julia
results = run_dust_optics(cfg)
df = annual_summary(results)
filter(:Varname => ==("ext2d_550nm_PM32Dust"), df)
```
"""
function run_dust_optics(cfg::DustOpticsConfig; runs=cfg.stats.runs)
    validate_optics_config(cfg)
    coefficients = optical_coefficients(cfg)
    return process_runs(r -> load_optics_run(cfg, coefficients, r), cfg.stats, runs)
end

"""
    optical_coefficients(cfg::DustOpticsConfig; write=true) -> Vector{OpticalCoefficients}

Optical properties of the size bins at the diagnostic wavelengths of every spectral set
of `cfg` ([`interpolate_coefficients`](@ref)); if `write`, they are also written to the
`coefficients_file` of each set that has one ([`write_coefficients`](@ref)).

# Example

```julia
coefs = optical_coefficients(cfg; write=false)
coefs[1].k_abs                             # (wavelength, bin)
```
"""
function optical_coefficients(cfg::DustOpticsConfig; write::Bool=true)
    coefs = OpticalCoefficients[]
    for set in cfg.spectral
        c = interpolate_coefficients(set, cfg.bins; extrapolation=cfg.extrapolation)
        if write && set.coefficients_file !== nothing
            write_coefficients(set.coefficients_file, c, cfg.bins)
            @info "Wrote the optical coefficients of $(set.name) to $(set.coefficients_file)"
        end
        push!(coefs, c)
    end
    return coefs
end
