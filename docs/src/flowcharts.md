# Flowcharts

Data flow through the two submodules, and the function calls within each. Rectangles are
functions, rounded boxes are data, cylinders are files. Functions from `PlanetTools` are
marked *(PlanetTools)*.

## 1. Data flow between the submodules

```mermaid
flowchart LR
    SA[("scaleacc monthly files<br/>&lt;MON&gt;&lt;YYYY&gt;.&lt;product&gt;&lt;run&gt;.nc<br/>aij, aijl, taijl, …")]
    TS["TimeSeriesMeSd<br/>run_timeseries_statistics"]
    O1[("time series and statistics<br/>e.g. Tracermass3d/…<br/>*.mts, *.ats, *.sts,<br/>*.moa/soa, *.mos/sos, *.mom/som")]
    MTS[("dust tracers q_s and air mass m_air<br/>MonthlyTs/*.ijl.mts.nc")]
    SP[("optical-property files aer_*.nc<br/>band bounds sp_*_bands.csv")]
    DO["DustOptics<br/>run_dust_optics"]
    O2[("time series and statistics<br/>e.g. Traceroptics3d/…<br/>r_eff, τ_sct, τ_abs, τ_ext")]
    CF[("interpolated coefficients<br/>aer_*_interpol_*.nc")]

    SA --> TS --> O1
    O1 -.- MTS
    MTS --> DO
    SP --> DO
    DO --> O2
    DO --> CF
```

## 2. Data flow within one run (shared pipeline)

`process_run(cfg, rd)` of `TimeSeriesMeSd`, used by both submodules. The variables are
processed one at a time; a `LazyVariable` is computed only here.

```mermaid
flowchart TD
    RD(["RunData<br/>grid, months, variables<br/>(VariableData or LazyVariable)"])
    OE["output_entries<br/>model levels; optionally pressure levels<br/>(interpolate_to_pressure, monthly p_s)"]
    X(["monthly field x<br/>(lon, lat[, level], time)"])
    SM["planet_means<br/>ij: gridded · zm: zonal mean · gm: global mean"]
    AM["annual_means<br/>full calendar years"]
    SE["seasonal_means<br/>one series per season"]
    MYA["multiyear_stats<br/>annual means in window"]
    MYS["multiyear_stats<br/>seasonal means in window"]
    MS["monthly_stats<br/>each month of the annual cycle"]
    WF["write_field!"]
    F1[("MonthlyTs/*.mts.nc")]
    F2[("AnnualTs/*.ats.nc")]
    F3[("SeasonalTs/*.sts.nc")]
    F4[("AnnualMeSd/*.moa.nc, *.soa.nc")]
    F5[("SeasonalMeSd/*.mos.nc, *.sos.nc")]
    F6[("MonthlyMeSd/*.mom.nc, *.som.nc")]
    RR(["RunResult<br/>files; annual and seasonal<br/>summaries of 2D global means"])

    RD --> OE --> X --> SM
    SM --> AM --> MYA
    SM --> SE --> MYS
    SM --> MS
    SM -- mts --> WF
    AM -- ats --> WF
    SE -- sts --> WF
    MYA -- moa, soa --> WF
    MYS -- mos, sos --> WF
    MS -- mom, som --> WF
    WF --> F1 & F2 & F3 & F4 & F5 & F6
    MYA --> RR
    MYS --> RR
```

## 3. Data flow within DustOptics

```mermaid
flowchart TD
    CFG(["DustOpticsConfig<br/>stats, bins, fractions, spectral sets"])
    SP[("aer_*.nc + sp_*_bands.csv")]
    RST["read_spectral_table"]
    ST(["SpectralTable → SpectralSet<br/>diagnostic wavelengths λ"])
    IC["interpolate_coefficients<br/>log–log in λ and r"]
    OC(["OpticalCoefficients<br/>k_scat, k_abs, g at (λ, r_b)"])
    CF[("aer_*_interpol_*.nc")]
    MTS[("tracer time series<br/>*.ijl.mts.nc: q_s, m_air")]
    LOR["load_optics_run"]
    DBM(["DustBinMasses<br/>m_b = f_b · q_s · m_air"])
    OV["optics_variables"]
    R3["effective_radius<br/>effRad3d"]
    RC["effective_radius ∘ column_masses<br/>effRadColMass"]
    RA["column_mean_effective_radius<br/>effRadColAve"]
    T3["combine_bins: τ = Σ k_b m_b<br/>sct3d, abs3d, ext3d"]
    T2["combine_bins ∘ column_masses<br/>sct2d, abs2d, ext2d"]
    PR["process_run<br/>(shared pipeline)"]

    SP --> RST --> ST --> IC --> OC
    CFG --> RST
    OC --> CF
    MTS --> LOR --> DBM
    CFG --> LOR
    DBM --> OV
    OC --> OV
    OV --> R3 & RC & RA & T3 & T2
    R3 & RC & RA & T3 & T2 --> PR
```

## 4. Function calls: TimeSeriesMeSd

```mermaid
flowchart LR
    RTS["run_timeseries_statistics"] --> VC["validate_config"]
    RTS --> PRs["process_run(cfg, run)"]
    PRs --> LR["load_run"]
    PRs --> PRd["process_run(cfg, rd)"]
    PRN["process_runs(load, cfg, runs)<br/>entry point for applications"] --> PRd

    LR --> INF["input_filename"]
    INF --> SCF["scaleacc_filename<br/>or cfg.filename"]
    SCF --> FN["fnames_modele<br/>(PlanetTools)"]
    LR --> SS["spatial_shape"]
    LR --> RC["read_coordinate"]

    PRd --> OE["output_entries"]
    OE --> OPL["output_pressure_levels"]
    OE --> ITP["interpolate_to_pressure"]
    ITP --> LP["layer_pressure"]
    OE --> FLD["field<br/>VariableData / LazyVariable"]

    PRd --> PER["annual_periods, seasonal_periods,<br/>climatology_slots, period_in_window"]
    PRd --> MW["month_weights"]
    PRd --> OF["output_filename"]
    PRd --> COF["create_output_file"]
    COF --> OPB["output_pressure_bounds"]
    OPB --> PLB["pressure_level_bounds"]
    OPB --> MLB["model_level_bounds"]
    COF --> WN["writenetcdf<br/>(PlanetTools)"]

    PRd --> PM["planet_means"]
    PRd --> AM["annual_means"]
    PRd --> SE["seasonal_means"]
    PRd --> MY["multiyear_stats"]
    PRd --> MS["monthly_stats"]
    PRd --> WF["write_field!"]
    WF --> WN

    PM --> CPM["calc_planet_means!<br/>(PlanetTools)"]
    AM --> CAM["calc_annual_mean_ts!<br/>(PlanetTools)"]
    SE --> CSM["calc_seasonal_mean_ts!<br/>(PlanetTools)"]
    MY --> CMS["calc_multiyear_stats!<br/>(PlanetTools)"]
    MS --> CMS

    PRd --> RRs["RunResult"]
    RRs --> SUM["annual_summary,<br/>seasonal_summary"]
```

## 5. Function calls: DustOptics

```mermaid
flowchart LR
    RDO["run_dust_optics"] --> VOC["validate_optics_config"]
    VOC --> OVN["optics_variable_names"]
    VOC --> VC["validate_config<br/>(TimeSeriesMeSd)"]
    RDO --> OC["optical_coefficients"]
    OC --> IC["interpolate_coefficients"]
    IC --> IEL["interpolate_extinction_coeff_log<br/>(PlanetTools)"]
    OC --> WC["write_coefficients"]
    WC --> WN["writenetcdf<br/>(PlanetTools)"]
    RDO --> PRN["process_runs<br/>(TimeSeriesMeSd)"]
    PRN --> LOR["load_optics_run"]
    PRN --> PRd["process_run(cfg, rd)<br/>(TimeSeriesMeSd)"]

    RSTb["read_spectral_table<br/>(application script)"] --> RBB["read_band_bounds"]
    RSTb --> ST(["SpectralTable → SpectralSet"])

    LOR --> IF["input_file"]
    LOR --> RF["read_field"]
    LOR --> RFR["run_fractions"]
    LOR --> BS["bin_sources"]
    LOR --> DBM(["DustBinMasses"])
    LOR --> OV["optics_variables"]
    OV --> OVS["optics_variable_specs"]
    OVS --> OCL["optics_classes"]
    OVS --> WL["wavelength_label"]
    OV --> LV(["LazyVariable per output variable"])

    LV -. "computed in process_run" .-> ER["effective_radius"]
    LV -.-> CMER["column_mean_effective_radius"]
    LV -.-> CB["combine_bins"]
    LV -.-> CM["column_masses"]
    ER --> CB
    CMER --> ER
    CMER --> CWM["calc_weighted_mean<br/>(PlanetTools)"]
```
