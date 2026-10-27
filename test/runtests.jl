using PlanetPostStats
using PlanetPostStats.TimeSeriesMeSd
using PlanetPostStats.TimeSeriesMeSd: absmonth, fromabsmonth, contiguous_range
using PlanetTools.PlanetStatistics: calc_planet_means!, calc_annual_mean_ts!
using Dates: Date, Year
using NCDatasets
using OrderedCollections: OrderedDict
using Statistics: mean, std
using Test

const NLON, NLAT, NLEV = 4, 3, 2
const LON = Float32[-135, -45, 45, 135]
const LAT = Float32[-60, 0, 60]
const AREA = [cosd(lat) for _ in LON, lat in LAT]

# Synthetic field values (physical units) as functions of cell, level, and month.
prsurf_val(i, j, ym) = 6.0 + ym.year + 0.01 * i
z_val(i, j, l, ym) = 1000.0 * l + j
airmass_val(i, j, l, ym) = Float64(ym.month)
clay_raw(i, j, l, ym) = 100.0 + ym.year          # file units: 1e-8 kg/kg
silt_raw(i, j, l, ym) = 50.0

"""Write a fake set of monthly scaleacc files for one run."""
function make_fake_run(
    root,
    run,
    months;
    skip=Dict{String,Vector{YearMonth}}(),
    missing_cells=Dict{YearMonth,Vector{Tuple{Int,Int}}}(),
)
    dir = joinpath(root, run, "Post")
    mkpath(dir)
    for ym in months
        for product in ("aij", "aijl", "taijl")
            ym in get(skip, product, YearMonth[]) && continue
            path = joinpath(dir, scaleacc_filename(ym, product, run))
            NCDataset(path, "c") do ds
                defVar(ds, "lon", LON, ("lon",))
                defVar(ds, "lat", LAT, ("lat",))
                defVar(
                    ds,
                    "axyp",
                    AREA,
                    ("lon", "lat");
                    attrib=OrderedDict("units" => "m2", "long_name" => "cell area"),
                )
                if product == "aij"
                    p = Float32[prsurf_val(i, j, ym) for i in 1:NLON, j in 1:NLAT]
                    for (i, j) in get(missing_cells, ym, Tuple{Int,Int}[])
                        p[i, j] = -1.0f30
                    end
                    defVar(
                        ds,
                        "prsurf",
                        p,
                        ("lon", "lat");
                        attrib=OrderedDict(
                            "units" => "hPa",
                            "long_name" => "surface pressure",
                            "missing_value" => -1.0f30,
                        ),
                    )
                else
                    defVar(ds, "level", Float32[1, 2], ("level",))
                    defVar(ds, "plm", Float32[5, 1], ("level",))
                    defVar(ds, "ple", Float32[3, 0.2], ("level",))
                    if product == "aijl"
                        defVar(
                            ds,
                            "z",
                            Float32[
                                z_val(i, j, l, ym) for i in 1:NLON, j in 1:NLAT, l in 1:NLEV
                            ],
                            ("lon", "lat", "level");
                            attrib=OrderedDict("units" => "m", "long_name" => "height"),
                        )
                        defVar(
                            ds,
                            "airmass",
                            Float32[
                                airmass_val(i, j, l, ym) for
                                i in 1:NLON, j in 1:NLAT, l in 1:NLEV
                            ],
                            ("lon", "lat", "level");
                            attrib=OrderedDict(
                                "units" => "kg/m2", "long_name" => "air mass"
                            ),
                        )
                    else
                        defVar(
                            ds,
                            "Clay",
                            Float32[
                                clay_raw(i, j, l, ym) for
                                i in 1:NLON, j in 1:NLAT, l in 1:NLEV
                            ],
                            ("lon", "lat", "level");
                            attrib=OrderedDict(
                                "units" => "1e-8 kg/kg", "long_name" => "Clay"
                            ),
                        )
                        defVar(
                            ds,
                            "Silt1",
                            Float32[
                                silt_raw(i, j, l, ym) for
                                i in 1:NLON, j in 1:NLAT, l in 1:NLEV
                            ],
                            ("lon", "lat", "level");
                            attrib=OrderedDict(
                                "units" => "1e-8 kg/kg", "long_name" => "Silt1"
                            ),
                        )
                    end
                end
            end
        end
    end
end

function test_config(root, runs; kw...)
    return TimeSeriesConfig(;
        simpath=root,
        runs=runs,
        groups=[
            InputGroup("aij", [InputVariable("prsurf"; units="hPa")]),
            InputGroup(
                "aijl",
                [
                    InputVariable("z"; units="m"),
                    InputVariable("airmass"; units="kg m-2", to_pressure=false),
                ],
            ),
            InputGroup(
                "taijl",
                [
                    InputVariable("Clay"; scale=1e-8, units="kg kg-1"),
                    InputVariable("Silt1"; scale=1e-8, units="kg kg-1"),
                ],
            ),
        ],
        derived=dust_mass_classes(["Clay", "Silt1"]),
        sim_period=(YearMonth(1, 1), YearMonth(4, 12)),
        stat_period=(YearMonth(2, 1), YearMonth(4, 12)),
        calendar=MARS_CALENDAR,
        seasons=MARS_SEASONS,
        vertical=VerticalCoordinate(31, 6.86, 0.2),
        output_dir="Out",
        output_name="test",
        global_attributes=OrderedDict{String,Any}("File_creator" => "test"),
        kw...,
    )
end

@testset "PlanetPostStats.TimeSeriesMeSd" begin
    @testset "calendar" begin
        @test YearMonth(1, 12) + 1 == YearMonth(2, 1)
        @test YearMonth(2, 1) - 1 == YearMonth(1, 12)
        @test fromabsmonth(absmonth(YearMonth(-3, 7))) == YearMonth(-3, 7)
        @test months_between(YearMonth(1, 1), YearMonth(20, 12)) == 239
        @test_throws ArgumentError YearMonth(1, 13)
        @test days_in_year(MARS_CALENDAR) == 669
        ms = yearmonth_range(YearMonth(1, 7), YearMonth(3, 12))
        @test [p.first.year for p in annual_periods(ms)] == [2, 3]
        ndj = seasonal_periods(ms, SeasonSpec("NDJ", 11, 3))
        @test [p.first for p in ndj] == [YearMonth(1, 11), YearMonth(2, 11)]
        @test season_months(SeasonSpec("DJF", 12, 3)) == [12, 1, 2]
        slots = climatology_slots(ms, (YearMonth(2, 1), YearMonth(3, 12)))
        @test slots[1] == [7, 19] && slots[12] == [18, 30]
        @test time_value(TimeAxis(), YearMonth(1, 1)) == 1
    end

    @testset "PlanetTools adapters" begin
        # hand-computed values (PlanetTools conventions)
        x = reshape(Float32[1, 3, 2, NaN], 2, 2, 1)
        sm = planet_means(x, [1.0 3.0; 1.0 3.0]; threshold=0)
        @test vec(sm.zm) == Float32[2, 2] && vec(sm.zm_w) == Float32[1, 0.5]
        @test sm.gm[1] ≈ 2 && sm.gm_w[1] ≈ 0.625
        sm80 = planet_means(x, [1.0 3.0; 1.0 3.0]; threshold=80)
        @test isnan(sm80.zm[1, 2]) && sm80.zm_w[1, 2] == 0
        me, sd, av = multiyear_stats(
            reshape(Float32[1, 2, 3, NaN], 1, 4), nothing, 1:4; threshold=50
        )
        @test me[1] ≈ 2 && sd[1] ≈ 1 && av[1] ≈ 0.75   # sample SD (n - 1)
        @test contiguous_range([3, 4, 5]) == 3:5
        @test_throws ArgumentError contiguous_range([3, 5])

        # the adapters (level chunking, reshapes) reproduce direct PlanetTools calls
        nlon, nlat, nlev, nt = 6, 4, 3, 36
        f = rand(Float32, nlon, nlat, nlev, nt)
        f[1:2, 2, 1, 5] .= NaN
        f[:, :, 2, 20] .= NaN
        area = [cosd(l) for _ in 1:nlon, l in (-60, -20, 20, 60)]
        ms = yearmonth_range(YearMonth(1, 1), YearMonth(3, 12))
        w = month_weights(MARS_CALENDAR, ms)
        sm = planet_means(f, area; threshold=80)
        F = reshape(f, 1, 1, size(f)...)
        zm = fill(NaN32, 1, 1, 1, nlat, nlev, nt)
        gm = fill(NaN32, 1, 1, 1, 1, nlev, nt)
        wz, wg = similar(zm), similar(gm)
        calc_planet_means!(F, zm, gm, wz, wg, area, [nt]; threshold=80)
        @test isequal(vec(zm), vec(sm.zm)) && isequal(vec(wg), vec(sm.gm_w))
        an = annual_means(f, sm, ms, w; threshold=80)
        A = fill(NaN32, 1, 1, nlon, nlat, nlev, 3)
        ZA = fill(NaN32, 1, 1, 1, nlat, nlev, 3)
        GA = fill(NaN32, 1, 1, 1, 1, nlev, 3)
        calc_annual_mean_ts!(
            F,
            zm,
            gm,
            A,
            ZA,
            GA,
            similar(A),
            similar(ZA),
            similar(GA),
            [Date(m.year, m.month) for m in ms],
            [nt],
            [Date(1):Year(1):Date(3)];
            zm_w=wz,
            gm_w=wg,
            month_weights=w,
            threshold=80,
        )
        @test isequal(vec(A), vec(an.ij)) && isequal(vec(ZA), vec(an.zm))
        @test isequal(vec(GA), vec(an.gm))
        se = seasonal_means(f, sm, ms, SeasonSpec("NDJ", 11, 3), nothing; threshold=80)
        @test size(se.gm) == (1, 1, nlev, 2)
        @test se.ij[1, 1, 1, 1] ≈ mean(f[1, 1, 1, 11:13])
        cm, cs, _ = monthly_stats(sm.gm, sm.gm_w, climatology_slots(ms, (ms[1], ms[end])))
        @test size(cm) == (1, 1, nlev, 12)
    end

    @testset "config validation" begin
        root = mktempdir()
        @test validate_config(test_config(root, ["A"])) isa TimeSeriesConfig
        @test_throws ArgumentError validate_config(test_config(root, ["A", "A"]))
        @test_throws ArgumentError validate_config(
            test_config(root, ["A"]; stat_period=(YearMonth(9, 1), YearMonth(9, 12)))
        )
        @test_throws ArgumentError validate_config(
            test_config(root, ["A"]; derived=[DerivedVariable("x", ["nope"], sum_of())])
        )
        @test_throws ArgumentError validate_config(test_config(root, ["A"]; threshold=0.0))
        pl = PressureLevels([6.0, 3.0])
        @test validate_config(test_config(root, ["A"]; pressure_levels=pl)) isa
            TimeSeriesConfig
        @test_throws ArgumentError validate_config(
            test_config(root, ["A"]; pressure_levels=pl, vertical=nothing)
        )
        @test_throws ArgumentError validate_config(
            test_config(
                root, ["A"]; pressure_levels=PressureLevels([6.0]; surface_pressure="ps")
            ),
        )
        @test_throws ArgumentError PressureLevels([6.0, -1.0])
        @test_throws ArgumentError PressureLevels([6.0, 6.0])
        @test_throws ArgumentError PressureLevels([6.0]; mode=:other)
    end

    @testset "pressure-level interpolation" begin
        vc = VerticalCoordinate(3, 10.0, 1.0)      # layers 1, 2 σ; layer 3 pressure
        @test layer_pressure(10.0, 9.0, 1, vc) ≈ 9.0
        @test layer_pressure(5.5, 9.0, 1, vc) ≈ 1.0 + 8.0 / 9.0 * 4.5
        @test layer_pressure(5.5, 0.5, 3, vc) == 0.5
        plm = [9.0, 5.0, 0.5]
        # x = ln p of the reference column: interpolation in ln p is exact
        x = reshape(Float32[log(9.0), log(5.0), log(0.5)], 1, 1, 3, 1)
        ps = fill(10.0f0, 1, 1, 1)
        y = interpolate_to_pressure(x, ps, plm, vc, [12.0, 9.5, 7.0, 2.0, 0.5, 0.1])
        @test size(y) == (1, 1, 6, 1)
        @test isnan(y[1])                       # below the surface
        @test y[2] ≈ log(9.0)                   # surface – lowest layer: lowest value
        @test y[3] ≈ log(7.0)
        @test y[4] ≈ log(2.0)
        @test y[5] ≈ log(0.5)
        @test isnan(y[6])                       # above the uppermost layer mid-point
        # lower surface pressure shifts the σ layers; NaN surface pressure → all missing
        ps2 = reshape(Float32[5.5, NaN], 2, 1, 1)
        x2 = repeat(x, 2)
        y2 = interpolate_to_pressure(x2, ps2, plm, vc, [7.0, 5.0])
        p1 = layer_pressure(5.5, 9.0, 1, vc)
        p2 = layer_pressure(5.5, 5.0, 2, vc)
        wgt = (log(5.0) - log(p1)) / (log(p2) - log(p1))
        @test isnan(y2[1, 1, 1, 1])
        @test y2[1, 1, 2, 1] ≈ log(9.0) + wgt * (log(5.0) - log(9.0)) rtol = 1e-6
        @test all(isnan, y2[2, 1, :, 1])
        @test_throws ArgumentError interpolate_to_pressure(x, ps, reverse(plm), vc, [5.0])
        @test_throws DimensionMismatch interpolate_to_pressure(
            x, fill(10.0f0, 2, 1, 1), plm, vc, [5.0]
        )
    end

    @testset "end-to-end" begin
        root = mktempdir()
        allmonths = yearmonth_range(YearMonth(1, 1), YearMonth(4, 12))
        # Run A: complete, except one missing taijl file and one missing aij cell.
        make_fake_run(
            root,
            "A",
            allmonths;
            skip=Dict("taijl" => [YearMonth(1, 5)]),
            missing_cells=Dict(YearMonth(2, 3) => [(1, 1)]),
        )
        # Run B: ends in month 6 of year 2 (shorter than the statistics window).
        make_fake_run(root, "B", yearmonth_range(YearMonth(1, 1), YearMonth(2, 6)))
        cfg = test_config(root, ["A", "B"])
        results = run_timeseries_statistics(cfg)
        resA, resB = results
        @test resA.months == allmonths
        @test last(resB.months) == YearMonth(2, 6)
        @test length(resA.files) == 54
        @test all(isfile, resA.files)

        # regression test for the original variable-slot shift bug: aijl variables must
        # be at their place even when the taijl file is missing
        rd = load_run(cfg, "A")
        t = findfirst(==(YearMonth(1, 5)), rd.months)
        @test rd.variables["z"].data[2, 3, 2, t] ≈ z_val(2, 3, 2, YearMonth(1, 5))
        @test all(isnan, rd.variables["Clay"].data[:, :, :, t])
        @test rd.variables["Clay"].data[1, 1, 1, 1] ≈ 101.0f-8
        @test rd.variables["PM2Dust"].data[1, 1, 1, 1] ≈ 101.0f-8
        @test rd.variables["PM4Dust"].data[1, 1, 1, 1] ≈ 151.0f-8
        @test !haskey(rd.variables, "PM8Dust")            # Silt2 not loaded
        @test rd.variables["PM4Dust"].attrib["units"] == "kg kg-1"
        @test isnan(rd.variables["prsurf"].data[1, 1, 15])     # missing_value decoded

        out = output_directory(cfg, "A")
        # multi-year annual means of the global mean surface pressure: years 2-4
        NCDataset(joinpath(out, "AnnualMeSd", "ANN_0002_0004.testA.gml.moa.nc")) do ds
            gm_prs = 6 + 3 + 0.01 * sum(AREA .* (1:NLON)) / sum(AREA)
            @test ds["prsurf"][1, 1] ≈ gm_prs rtol = 1e-4
            @test ds["LS1"][] == 31
            @test ds["Ptop"][] ≈ 0.2f0
            @test ds["Clay"].attrib["units"] == "kg kg-1"
            @test !haskey(ds["Clay"].attrib, "scale_factor")
            @test ds["prsurf"].attrib["cell_methods"] ==
                "time: multi-year annual mean; lon, lat: mean"
        end
        NCDataset(joinpath(out, "AnnualMeSd", "ANN_0002_0004.testA.ijl.soa.nc")) do ds
            @test ds["prsurf"][2, 2] ≈ 1.0 rtol = 1e-5   # sample SD of 3 years
        end
        NCDataset(joinpath(out, "AnnualMeSd", "ANN_0002_0004.testA.gml.soa.nc")) do ds
            @test ds["Clay"][1, 1, 1] ≈ 1.0f-8 rtol = 1e-4
        end
        # annual mean of airmass (= month number) is the day-weighted mean month
        w = MARS_CALENDAR.days_in_month ./ 669
        NCDataset(joinpath(out, "AnnualTs", "ANN_0001-0004.testA.ijl.ats.nc")) do ds
            @test ds["airmass"][1, 1, 1, 1] ≈ sum(w .* (1:12)) rtol = 1e-6
            @test NCDatasets.variable(ds, "time")[:] == [1, 13, 25, 37]
            @test NCDatasets.variable(ds, "time").attrib["calendar"] == "mars"
            @test NCDatasets.variable(ds, "time").attrib["month_lengths"] ==
                collect(MARS_CALENDAR.days_in_month)
            @test ds["year"][:] == 1:4
            @test ds["Dates"][:] == ["0001-01-01", "0002-01-01", "0003-01-01", "0004-01-01"]
        end
        # seasonal mean FMA, weighted by month length
        d = MARS_CALENDAR.days_in_month
        fma = (d[2] * 2 + d[3] * 3 + d[4] * 4) / (d[2] + d[3] + d[4])
        NCDataset(joinpath(out, "SeasonalTs", "FMA_0001-0004.testA.zml.sts.nc")) do ds
            @test ds["airmass"][1, 2, 1, 1] ≈ fma rtol = 1e-6
            @test size(ds["airmass"]) == (1, NLAT, NLEV, 4)
        end
        # NDJ season wraps into the next year: only 3 complete occurrences in years 1-4
        @test isfile(joinpath(out, "SeasonalTs", "NDJ_0001-0003.testA.ijl.sts.nc"))
        # monthly climatology
        NCDataset(joinpath(out, "MonthlyMeSd", "JAN-DEC_0002_0004.testA.gml.mom.nc")) do ds
            @test ds["airmass"][1, 1, 1, :] ≈ 1:12
        end
        NCDataset(joinpath(out, "MonthlyMeSd", "JAN-DEC_0002_0004.testA.gml.som.nc")) do ds
            @test all(iszero, ds["airmass"][1, 1, 1, :])
        end
        # fill values are written as _FillValue and read back as missing
        NCDataset(joinpath(out, "MonthlyTs", "JAN0001-DEC0004.testA.ijl.mts.nc")) do ds
            @test ismissing(ds["prsurf"][1, 1, 15])
            @test ismissing(ds["Clay"][1, 1, 1, 5])
            @test ds["prsurf"].attrib["_FillValue"] == -1.0f30
        end

        # summary tables
        ann = annual_summary(results)
        @test only(ann[ann.Experiment .== "A", :Varname]) == "prsurf"
        @test isempty(ann[ann.Experiment .== "B", :])   # no full year 2 in B
        sea = seasonal_summary(results)
        @test count(==("A"), sea.Experiment) == 4

        # unweighted seasons reproduce the original arithmetic seasonal means
        cfgu = test_config(
            root, ["A"]; weight_seasons=false, products=[:sts], output_dir="OutU"
        )
        run_timeseries_statistics(cfgu)
        NCDataset(
            joinpath(
                output_directory(cfgu, "A"), "SeasonalTs", "FMA_0001-0004.testA.gml.sts.nc"
            ),
        ) do ds
            @test ds["airmass"][1, 1, 1, 1] ≈ 3.0
        end
    end

    @testset "pressure levels end-to-end" begin
        root = mktempdir()
        allmonths = yearmonth_range(YearMonth(1, 1), YearMonth(4, 12))
        make_fake_run(root, "A", allmonths)
        plev = [8.0, 6.0, 3.0, 1.5, 0.5]
        vc = VerticalCoordinate(31, 6.86, 0.2)
        plm = Float32[5, 1]
        # expected value of z in column (i, j), month ym, at pressure P (both layers σ)
        function z_expected(i, j, ym, P)
            ps = Float32(prsurf_val(i, j, ym))
            p1 = layer_pressure(ps, plm[1], 1, vc)
            p2 = layer_pressure(ps, plm[2], 2, vc)
            (P > ps || P < p2) && return NaN
            P >= p1 && return z_val(i, j, 1, ym)
            wgt = (log(P) - log(p1)) / (log(p2) - log(p1))
            return z_val(i, j, 1, ym) + wgt * (z_val(i, j, 2, ym) - z_val(i, j, 1, ym))
        end

        # mode :add — additional files with tags ijp, zmp, gmp
        cfg = test_config(
            root,
            ["A"];
            pressure_levels=PressureLevels(plev),
            products=[:mts, :ats, :moa],
            output_dir="OutAdd",
        )
        res = only(run_timeseries_statistics(cfg))
        @test length(res.files) == 2 * 3 * 4
        @test all(isfile, res.files)
        out = output_directory(cfg, "A")
        NCDataset(joinpath(out, "MonthlyTs", "JAN0001-DEC0004.testA.ijp.mts.nc")) do ds
            @test ds["plev"][:] ≈ plev
            @test ds["plev"].attrib["units"] == "hPa"
            @test ds["plev"].attrib["positive"] == "down"
            @test haskey(ds.attrib, "Vertical_interpolation")
            @test !haskey(ds, "airmass")         # to_pressure = false
            @test !haskey(ds, "prsurf")          # 2D
            @test !haskey(ds, "plm")
            @test dimnames(ds["z"]) == ("lon", "lat", "plev", "time")
            for (t, ym) in enumerate(allmonths),
                k in eachindex(plev),
                (i, j) in ((1, 1), (4, 3))

                zexp = z_expected(i, j, ym, plev[k])
                if isnan(zexp)
                    @test ismissing(ds["z"][i, j, k, t])
                else
                    @test ds["z"][i, j, k, t] ≈ zexp rtol = 1e-5
                end
            end
            # constant-in-height mixing ratios stay constant
            @test ds["Clay"][2, 2, 3, 1] ≈ 101.0f-8 rtol = 1e-5
            @test ds["PM4Dust"][2, 2, 3, 1] ≈ 151.0f-8 rtol = 1e-5
        end
        # model-level files unchanged
        NCDataset(joinpath(out, "MonthlyTs", "JAN0001-DEC0004.testA.ijl.mts.nc")) do ds
            @test !haskey(ds, "plev")
            @test haskey(ds, "airmass")
            @test ds["z"][2, 3, 2, 1] ≈ z_val(2, 3, 2, YearMonth(1, 1))
        end
        # 8 hPa is below the surface in year 1 (ps ≈ 7 hPa), above it in years 2–4
        NCDataset(joinpath(out, "AnnualTs", "ANN_0001-0004.testA.gmp.ats.nc")) do ds
            @test ismissing(ds["z"][1, 1, 1, 1])
            @test !ismissing(ds["z"][1, 1, 1, 2])
            @test ds["z"][1, 1, 3, 1] > ds["z"][1, 1, 2, 1]
        end

        # mode :replace — pressure-level fields in the regular files
        cfgr = test_config(
            root,
            ["A"];
            pressure_levels=PressureLevels(plev; mode=:replace),
            products=[:mts],
            output_dir="OutRep",
        )
        resr = only(run_timeseries_statistics(cfgr))
        @test length(resr.files) == 3
        NCDataset(
            joinpath(
                output_directory(cfgr, "A"), "MonthlyTs", "JAN0001-DEC0004.testA.zml.mts.nc"
            ),
        ) do ds
            @test dimnames(ds["z"]) == ("lon", "lat", "plev", "time")
            @test dimnames(ds["airmass"]) == ("lon", "lat", "level", "time")
            @test dimnames(ds["prsurf"]) == ("lon", "lat", "time")
            @test haskey(ds, "plm")
        end
    end
end
