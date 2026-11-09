# Tests of PlanetPostStats.DustOptics; included from runtests.jl (uses its fixture).

using PlanetPostStats.DustOptics
using PlanetTools.PlanetAerosolUtils: calc_size_integrated_effective_radius

# Optical table following power laws, k = A λ^a r^b: interpolation linear in log-log
# space is exact, also when extrapolating linearly.
ksca_law(λ, r) = 2.0e3 * (λ / 1e-6)^-0.5 * (r / 1e-6)^-1.0
kabs_law(λ, r) = 3.0e2 * (λ / 1e-6)^0.7 * (r / 1e-6)^-0.8
gasy_law(λ, r) = 0.7 * (λ / 1e-6)^-0.1 * (r / 1e-6)^0.05

function power_law_table()
    bounds = [(0.3e-6, 0.5e-6), (0.5e-6, 0.8e-6), (0.8e-6, 1.2e-6), (1.2e-6, 2.0e-6)]
    centers = [sqrt(a * b) for (a, b) in bounds]
    radii = [0.1e-6, 0.3e-6, 1.0e-6, 3.0e-6]
    tab(f) = [f(λ, r) for λ in centers, r in radii]
    return SpectralTable(bounds, radii, tab(ksca_law), tab(kabs_law), tab(gasy_law))
end

const OPT_BINS = [
    DustBin("Clay1"; source="Clay", radius=0.132e-6, density=2.5e3),
    DustBin("Clay2"; source="Clay", radius=0.23e-6, density=2.5e3),
    DustBin("Clay3"; source="Clay", radius=0.416e-6, density=2.5e3),
    DustBin("Clay4"; source="Clay", radius=0.766e-6, density=2.5e3),
    DustBin("Silt1"; radius=1.386e-6, density=2.65e3),
]
const CLAY_SPLIT = [0.009, 0.081, 0.234, 0.676]
const OPT_WAVELENGTHS = [0.55e-6, 0.89e-6]

@testset "PlanetPostStats.DustOptics" begin
    @testset "bins and kernels" begin
        @test bin_sources(OPT_BINS) == ["Clay", "Silt1"]
        f = bin_fractions(OPT_BINS, "Clay" => CLAY_SPLIT)
        @test f == [CLAY_SPLIT; 1.0]
        @test_throws ArgumentError bin_fractions(OPT_BINS, "Clay" => [0.5, 0.5])
        @test_throws ArgumentError bin_fractions(OPT_BINS, "Sand" => [1.0])
        @test_throws ArgumentError DustBin("x"; radius=0, density=2500)
        @test optics_classes(OPT_BINS) == [
            "SumClay1Dust" => 1,
            "SumClay2Dust" => 2,
            "SumClay3Dust" => 3,
            "PM2Dust" => 4,
            "PM4Dust" => 5,
        ]
        @test isempty(optics_classes([DustBin("Silt1"; radius=1e-6, density=2650)]))

        sz = (3, 2, 4, 5)
        clay = rand(Float32, sz...) .* 1.0f-3
        silt = rand(Float32, sz...) .* 1.0f-3
        clay[1, 1, 1, 1] = NaN32
        clay[2, 1, :, 2] .= 0
        silt[2, 1, :, 2] .= 0
        m = DustBinMasses(["Clay", "Silt1"], [clay, silt], [1, 1, 1, 1, 2], f)
        # combine_bins against the explicit sum over the bins
        c = [1.0, 2.0, 3.0, 4.0, 5.0]
        expl = sum(c[b] * f[b] .* (b <= 4 ? clay : silt) for b in 1:5)
        @test combine_bins(m, c) ≈ expl nans = true
        # a zero coefficient does not propagate missing values of its source
        @test !isnan(combine_bins(m, [0, 0, 0, 0, 1.0])[1, 1, 1, 1])
        @test isnan(combine_bins(m, [1.0, 0, 0, 0, 1.0])[1, 1, 1, 1])
        # effective radius against PlanetTools, cell by cell
        reff = effective_radius(m, OPT_BINS)
        rho = [b.density for b in OPT_BINS]
        rad = [b.radius for b in OPT_BINS]
        for I in CartesianIndices(sz)
            mb = [f[b] * (b <= 4 ? clay[I] : silt[I]) for b in 1:5]
            if any(isnan, mb) || sum(mb) == 0
                @test isnan(reff[I])
            else
                @test reff[I] ≈ calc_size_integrated_effective_radius(mb, rho, rad) rtol =
                    1e-5
            end
        end
        # density-scaled opacity: τ / m_air with the masses m = f q m_air
        air = 1.0f0 .+ rand(Float32, sz...)
        air[3, 2, 1, 1] = 0
        mq = DustBinMasses(
            ["Clay", "Silt1"], [clay .* air, silt .* air], [1, 1, 1, 1, 2], f
        )
        dso = density_scaled_opacity(mq, c, air)
        expq = sum(c[b] * f[b] .* (b <= 4 ? clay : silt) for b in 1:5)
        @test dso[2, 2, 3, 4] ≈ expq[2, 2, 3, 4] rtol = 1e-5
        @test isnan(dso[1, 1, 1, 1]) && isnan(dso[3, 2, 1, 1])
        @test_throws DimensionMismatch density_scaled_opacity(mq, c, air[:, :, 1:2, :])
        # columns
        cm = column_masses(m)
        @test size(cm.masses[1]) == (3, 2, 1, 5)
        @test cm.masses[2][:, :, 1, :] ≈ dropdims(sum(silt; dims=3); dims=3)
        @test isnan(cm.masses[1][1, 1, 1, 1])
        # mass-weighted column average of the effective radius
        rbar = column_mean_effective_radius(m, OPT_BINS)
        @test size(rbar) == (3, 2, 5)
        mtot = combine_bins(m, ones(5))
        i, j, t = 2, 2, 3
        @test rbar[i, j, t] ≈
            sum(mtot[i, j, :, t] .* reff[i, j, :, t]) / sum(mtot[i, j, :, t]) rtol = 1e-5
        @test isnan(rbar[2, 1, 2])                    # no dust in the column
        @test !isnan(rbar[1, 1, 1])                   # missing layer ignored
    end

    @testset "spectral tables" begin
        t = power_law_table()
        set = SpectralSet("sw", t, OPT_WAVELENGTHS)
        @test wavelength_label(0.55e-6) == "550nm"
        @test wavelength_label(21.6e-6) == "21600nm"
        @test_throws ArgumentError SpectralSet("sw", t, [0.55e-6, 0.5500001e-6])
        c = interpolate_coefficients(set, OPT_BINS)
        for (w, λ) in enumerate(OPT_WAVELENGTHS), (b, bin) in enumerate(OPT_BINS)
            @test c.k_scat[w, b] ≈ ksca_law(λ, bin.radius) rtol = 1e-10
            @test c.k_abs[w, b] ≈ kabs_law(λ, bin.radius) rtol = 1e-10
            @test c.g_asym[w, b] ≈ gasy_law(λ, bin.radius) rtol = 1e-10
        end

        # files: band bounds (CSV, header in line 3, bands in reverse order of their
        # wavelengths, bounds in either order) and optical properties (NetCDF)
        dir = mktempdir()
        csv = joinpath(dir, "bands.csv")
        open(csv, "w") do io
            println(io, "# spectral bands of a test file")
            println(io, "# wavelengths in m")
            println(io, "Band,Lower,Upper")
            for (n, k) in enumerate(reverse(eachindex(t.band_bounds)))
                lo, hi = t.band_bounds[k]
                println(io, n, ",", n == 1 ? hi : lo, ",", n == 1 ? lo : hi)
            end
        end
        @test read_band_bounds(csv) == reverse(t.band_bounds)
        nc = joinpath(dir, "aer.nc")
        nb, nr = length(t.band_bounds), length(t.radius)
        rev = reverse(1:nb)
        NCDataset(nc, "c") do ds
            defVar(ds, "radius_eff_dry", repeat(t.radius, 1, 2), ("radius", "component"))
            for (v, x) in (("k_scat", t.k_scat), ("k_abs", t.k_abs), ("g_asym", t.g_asym))
                a = zeros(nb, 2, nr, 2)
                a[:, 1, :, 1] .= x[rev, :]
                a[:, 2, :, 1] .= 2 .* x[rev, :]          # other humidity
                defVar(ds, v, a, ("band", "humidity", "radius", "component"))
            end
        end
        tf = read_spectral_table(nc, csv)
        @test tf.band_bounds == reverse(t.band_bounds)
        @test tf.k_scat ≈ t.k_scat[rev, :]
        @test tf.radius ≈ t.radius
        cf = interpolate_coefficients(SpectralSet("sw", tf, OPT_WAVELENGTHS), OPT_BINS)
        @test cf.k_abs ≈ c.k_abs rtol = 1e-10            # band order does not matter
        @test read_spectral_table(nc, csv; humidity_index=2).k_abs ≈ 2 .* tf.k_abs
        bad = joinpath(dir, "bad.csv")
        write(bad, "a\nb\nBand,Lower,Upper\n1,1e-7,2e-7\n")
        @test_throws DimensionMismatch read_spectral_table(nc, bad)

        out = write_coefficients(joinpath(dir, "sub", "coef.nc"), c, OPT_BINS)
        NCDataset(out) do ds
            @test ds["wavelength"][:] ≈ OPT_WAVELENGTHS
            @test ds["radius_eff_dry"][:] ≈ [b.radius for b in OPT_BINS]
            @test ds["diag_k_scat_sw"][:, :] ≈ c.k_scat rtol = 1e-6
            @test ds["diag_g_asym_sw"].attrib["long_name"] == "Asymmetry parameter"
        end
    end

    @testset "end-to-end" begin
        root = mktempdir()
        allmonths = yearmonth_range(YearMonth(1, 1), YearMonth(4, 12))
        make_fake_run(root, "A", allmonths; skip=Dict("taijl" => [YearMonth(1, 5)]))
        # dust tracer time series with TimeSeriesMeSd
        tracer = test_config(
            root,
            ["A"];
            products=[:mts],
            output_dir="Tracermass3d",
            output_name="tracermass3d",
            derived=DerivedVariable[],
        )
        run_timeseries_statistics(tracer)

        coefdir = mktempdir()
        set = SpectralSet(
            "sw",
            power_law_table(),
            OPT_WAVELENGTHS;
            coefficients_file=joinpath(coefdir, "coef_sw.nc"),
        )
        stats = test_config(
            root,
            ["A"];
            groups=InputGroup[],
            derived=DerivedVariable[],
            products=[:mts, :ats, :moa],
            output_dir="Traceroptics3d",
            output_name="traceroptics3d",
        )
        cfg = DustOpticsConfig(;
            stats=stats,
            bins=OPT_BINS,
            spectral=[set],
            fractions=Dict("A" => bin_fractions(OPT_BINS, "Clay" => CLAY_SPLIT)),
            input_dir="Tracermass3d",
            input_name="tracermass3d",
        )
        names = optics_variable_names(cfg)
        # 3 effective radii + 3 quantities × 2 wavelengths × 3 (layers, columns,
        # density-scaled opacities) × (5 bins + 5 classes)
        @test length(names) == 3 + 3 * 2 * 3 * 10
        @test "extdso_550nm_Clay1" in names && "absdso_890nm_PM4Dust" in names
        cfgn = DustOpticsConfig(;
            stats=stats,
            bins=OPT_BINS,
            spectral=[set],
            density_scaled=false,
            input_dir="Tracermass3d",
            input_name="tracermass3d",
        )
        @test length(optics_variable_names(cfgn)) == 3 + 3 * 2 * 2 * 10
        @test !any(contains("dso_"), optics_variable_names(cfgn))
        m1 = DustBinMasses(
            ["Clay", "Silt1"],
            [ones(Float32, 1, 1, 1, 1), ones(Float32, 1, 1, 1, 1)],
            [1, 1, 1, 1, 2],
            ones(5),
        )
        @test_throws ArgumentError optics_variables(
            cfg, optical_coefficients(cfg; write=false), m1
        )
        @test length(optics_variables(cfgn, optical_coefficients(cfgn; write=false), m1)) ==
            3 + 3 * 2 * 2 * 10
        @test names[1:3] ==
            ["effRad3d_PM4Dust", "effRadColMass_PM4Dust", "effRadColAve_PM4Dust"]
        @test "ext3d_550nm_Clay1" in names && "sct2d_890nm_PM4Dust" in names
        @test validate_optics_config(cfg) === cfg
        @test endswith(input_file(cfg, "A"), "JAN0001-DEC0004.tracermass3dA.ijl.mts.nc")

        results = run_dust_optics(cfg)
        res = only(results)
        @test res.months == allmonths
        @test all(isfile, res.files)
        @test length(res.files) == 3 * 4      # (mts, ats, moa, soa) × (ij, zm, gm)
        @test isfile(joinpath(coefdir, "coef_sw.nc"))

        # expected values: mass of bin b = f_b q air, q of the source, air = month
        f = [CLAY_SPLIT; 1.0]
        q(b, ym) = b <= 4 ? clay_raw(1, 1, 1, ym) * 1e-8 : silt_raw(1, 1, 1, ym) * 1e-8
        mb(b, ym) = f[b] * q(b, ym) * airmass_val(1, 1, 1, ym)
        kext(w, b) =
            ksca_law(OPT_WAVELENGTHS[w], OPT_BINS[b].radius) +
            kabs_law(OPT_WAVELENGTHS[w], OPT_BINS[b].radius)
        rho = [b.density for b in OPT_BINS]
        rad = [b.radius for b in OPT_BINS]
        reff(ym) = calc_size_integrated_effective_radius([mb(b, ym) for b in 1:5], rho, rad)

        out = output_directory(stats, "A")
        ym = YearMonth(2, 7)
        t = findfirst(==(ym), allmonths)
        NCDataset(
            joinpath(out, "MonthlyTs", "JAN0001-DEC0004.traceroptics3dA.ijl.mts.nc")
        ) do ds
            @test ds["ext3d_550nm_Clay1"][2, 3, 1, t] ≈ kext(1, 1) * mb(1, ym) rtol = 1e-5
            @test ds["ext3d_890nm_PM4Dust"][1, 2, 2, t] ≈
                sum(kext(2, b) * mb(b, ym) for b in 1:5) rtol = 1e-5
            @test ds["sct3d_550nm_SumClay2Dust"][1, 1, 1, t] ≈
                sum(ksca_law(0.55e-6, rad[b]) * mb(b, ym) for b in 1:2) rtol = 1e-5
            @test ds["abs2d_550nm_Silt1"][3, 1, t] ≈
                NLEV * kabs_law(0.55e-6, rad[5]) * mb(5, ym) rtol = 1e-5
            @test ds["ext2d_550nm_PM4Dust"][4, 3, t] ≈
                NLEV * sum(kext(1, b) * mb(b, ym) for b in 1:5) rtol = 1e-5
            @test ds["effRad3d_PM4Dust"][1, 1, 1, t] ≈ reff(ym) rtol = 1e-5
            # density-scaled opacity: k f q, independent of the air mass
            @test ds["extdso_550nm_Clay1"][2, 3, 2, t] ≈
                kext(1, 1) * mb(1, ym) / airmass_val(1, 1, 1, ym) rtol = 1e-5
            @test ds["sctdso_890nm_PM4Dust"][1, 1, 1, t] ≈
                sum(ksca_law(0.89e-6, rad[b]) * f[b] * q(b, ym) for b in 1:5) rtol = 1e-5
            @test ds["extdso_550nm_Clay1"].attrib["units"] == "m2 kg-1"
            @test ds["extdso_550nm_PM4Dust"].attrib["long_name"] ==
                "PM4Dust Extinction Density-Scaled Opacity 550 nm"
            @test dimnames(ds["extdso_550nm_PM4Dust"]) == ("lon", "lat", "level", "time")
            @test ds["effRadColMass_PM4Dust"][1, 1, t] ≈ reff(ym) rtol = 1e-5
            @test ds["effRadColAve_PM4Dust"][1, 1, t] ≈ reff(ym) rtol = 1e-5
            @test ds["effRad3d_PM4Dust"].attrib["units"] == "m"
            @test ds["ext2d_550nm_PM4Dust"].attrib["long_name"] ==
                "PM4Dust Extinction AOT 550 nm"
            @test ds["ext3d_550nm_Clay1"].attrib["long_name"] ==
                "Extinction AOD 550 nm Clay1"
            @test dimnames(ds["ext2d_550nm_PM4Dust"]) == ("lon", "lat", "time")
            # month with the missing dust file: missing in the output
            t5 = findfirst(==(YearMonth(1, 5)), allmonths)
            @test ismissing(ds["ext3d_550nm_Clay1"][1, 1, 1, t5])
            @test ismissing(ds["effRadColAve_PM4Dust"][1, 1, t5])
            @test !haskey(ds, "Clay") && !haskey(ds, "airmass")
            @test haskey(ds, "plm") && haskey(ds, "LS1")
        end
        ann = annual_summary(results)
        @test "ext2d_550nm_PM4Dust" in ann.Varname
        @test "effRadColMass_PM4Dust" in ann.Varname
        @test !("ext3d_550nm_Clay1" in ann.Varname)

        # pressure levels: only the effective radius of the layers is interpolated
        statsp = test_config(
            root,
            ["A"];
            groups=InputGroup[],
            derived=DerivedVariable[],
            products=[:mts],
            output_dir="TraceropticsP",
            output_name="traceroptics3d",
            pressure_levels=PressureLevels([6.0, 3.0]),
        )
        cfgp = DustOpticsConfig(;
            stats=statsp,
            bins=OPT_BINS,
            spectral=[set],
            fractions=run -> bin_fractions(OPT_BINS, "Clay" => CLAY_SPLIT),
            quantities=[:ext],
            input_dir="Tracermass3d",
            input_name="tracermass3d",
        )
        resp = only(run_dust_optics(cfgp))
        @test length(resp.files) == 6
        NCDataset(
            joinpath(
                output_directory(statsp, "A"),
                "MonthlyTs",
                "JAN0001-DEC0004.traceroptics3dA.ijp.mts.nc",
            ),
        ) do ds
            @test haskey(ds, "effRad3d_PM4Dust")
            @test dimnames(ds["effRad3d_PM4Dust"]) == ("lon", "lat", "plev", "time")
            @test !haskey(ds, "ext3d_550nm_Clay1") && !haskey(ds, "prsurf")
            # density-scaled opacities are interpolated to the pressure levels
            @test dimnames(ds["extdso_550nm_PM4Dust"]) == ("lon", "lat", "plev", "time")
            @test ds["extdso_550nm_PM4Dust"][1, 1, 2, t] ≈
                sum(kext(1, b) * f[b] * q(b, ym) for b in 1:5) rtol = 1e-5
            @test !haskey(ds, "ext2d_550nm_PM4Dust")
            # vertically constant composition: the interpolated radius is unchanged
            @test ds["effRad3d_PM4Dust"][1, 1, 2, t] ≈ reff(ym) rtol = 1e-5
        end

        # input file selection
        @test_throws ArgumentError input_file(
            DustOpticsConfig(;
                stats=stats, bins=OPT_BINS, input_dir="Nope", input_name="tracermass3d"
            ),
            "A",
        )
        mts = dirname(input_file(cfg, "A"))
        cp(input_file(cfg, "A"), joinpath(mts, "JAN0001-DEC0003.tracermass3dA.ijl.mts.nc"))
        @test_throws ArgumentError input_file(cfg, "A")
        cfgx = DustOpticsConfig(;
            stats=stats,
            bins=OPT_BINS,
            input_dir="Tracermass3d",
            input_name="tracermass3d",
            input_prefix="JAN0001-DEC0004",
        )
        @test isfile(input_file(cfgx, "A"))

        # configuration checks
        @test_throws ArgumentError validate_optics_config(
            DustOpticsConfig(;
                stats=stats,
                bins=OPT_BINS,
                quantities=[:foo],
                input_dir="Tracermass3d",
                input_name="tracermass3d",
            ),
        )
        # without fractions: every bin is the whole source tracer
        @test run_fractions(
            DustOpticsConfig(;
                stats=stats,
                bins=[DustBin("Clay"; radius=1e-6, density=2500)],
                input_dir="Tracermass3d",
                input_name="tracermass3d",
            ),
            "A",
        ) == [1.0]
        @test_throws ArgumentError run_fractions(
            DustOpticsConfig(;
                stats=stats,
                bins=OPT_BINS,
                fractions=[1.0],
                input_dir="Tracermass3d",
                input_name="tracermass3d",
            ),
            "A",
        )
        @test_logs (:warn, r"sum to") run_fractions(
            DustOpticsConfig(;
                stats=stats,
                bins=OPT_BINS,
                fractions=[0.5, 0.1, 0.1, 0.1, 1.0],
                input_dir="Tracermass3d",
                input_name="tracermass3d",
            ),
            "A",
        )
    end
end
