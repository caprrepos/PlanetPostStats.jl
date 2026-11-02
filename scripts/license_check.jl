## Author: Jan P. Perlwitz, Climate, Aerosol, and Pollution Research, LLC
##         (jan.p.perlwitz@caprllc.com)
## AI cowork: Claude Opus 5.5 (model ID claude-opus-5-5, Anthropic),
##   used from 2026-10-24 to 2026-10-31

# License inventory of all dependencies of PlanetPostStats.
#
# Resolves the package with its weak (extension) and test dependencies in a temporary
# environment and reports them in separate scopes. It finds the license file of every
# Julia package and of every binary artifact of the `_jll` packages for this platform, and
# identifies the licenses with LicenseCheck.jl (SPDX license list).
# Writes a Markdown table and exits with status 1 if a license is copyleft or could not be
# identified, so the script can also run in CI.
#
# Usage, from the repository root:
#     julia scripts/license_check.jl [THIRD_PARTY_LICENSES.md]
# The binary artifacts are those of the host platform. To audit another platform, set
# LICENSE_CHECK_PLATFORM to its triplet; its artifacts are then downloaded, e.g.
#     LICENSE_CHECK_PLATFORM=aarch64-apple-darwin julia scripts/license_check.jl out.md

# No precompilation is needed to read license files. Disabling Pkg's automatic
# precompilation also avoids the "Unknown package" failure of `Pkg.add` in Julia 1.13.2
# (Pkg.jl issue #4866).
ENV["JULIA_PKG_PRECOMPILE_AUTO"] = "0"

using Pkg

const ROOT = dirname(@__DIR__)
const OUTFILE = isempty(ARGS) ? joinpath(ROOT, "THIRD_PARTY_LICENSES.md") : ARGS[1]
const PROJECT = Pkg.TOML.parsefile(joinpath(ROOT, "Project.toml"))

Pkg.activate(; temp=true, io=devnull)
Pkg.develop(; path=ROOT, io=devnull)
const WEAKDEPS = collect(keys(get(PROJECT, "weakdeps", Dict())))
const EXTRAS = collect(keys(get(PROJECT, "extras", Dict())))
Pkg.add(unique([WEAKDEPS; EXTRAS; "LicenseCheck"]); io=devnull)
Pkg.instantiate(; io=devnull)

using Artifacts: artifact_meta
using Base.BinaryPlatforms: HostPlatform, triplet
using LicenseCheck: find_licenses, licensecheck

const PERMISSIVE = Set([
    "MIT",
    "MIT-0",
    "X11",
    "ISC",
    "BSD-2-Clause",
    "BSD-3-Clause",
    "BSD-4-Clause-UC",
    "Zlib",
    "Apache-2.0",
    "BSL-1.0",
    "curl",
    "Unlicense",
    "CC0-1.0",
    "libpng-2.0",
    "Libpng",
    "PSF-2.0",
    "Python-2.0",
    "OpenSSL",
    "SSLeay-standalone",
    "IJG",
    "FTL",
    "NCSA",
    "BSD-Source-Code",
    "HPND",
    "MIT-CMU",
    "bzip2-1.0.6",
    "LLVM-exception",
    "0BSD",
    "BSD-3-Clause-Clear",
    "BSD-3-Clause-Open-MPI",
    "Info-ZIP",
    "X11-distribute-modifications-variant",
])
const WEAK_COPYLEFT = r"^(LGPL|MPL|EPL|CDDL|GPL-2\.0-with-classpath)"
const STRONG_COPYLEFT = r"^(GPL|AGPL)"

# Manually reviewed entries (2026-10-24): class and finding. Entries not listed here that
# are not clearly permissive make the script fail, so that new dependencies are reviewed.
const REVIEWED = Dict(
    "FFMPEG_jll → FFMPEG binary" => (
        "strong copyleft, nonfree",
        "built with --enable-gpl --enable-nonfree --enable-libfdk-aac; `ffmpeg -L`: " *
        "\"not legally redistributable\"; via Makie",
    ),
    "x264_jll → x264 binary" => ("strong copyleft", "GPL-2.0-or-later; via Makie → FFMPEG"),
    "x265_jll → x265 binary" => ("strong copyleft", "GPL-2.0-or-later; via Makie → FFMPEG"),
    "libfdk_aac_jll → libfdk_aac binary" => (
        "custom",
        "Fraunhofer FDK AAC license; FFmpeg treats it as nonfree with GPL; via Makie",
    ),
    "Rmath_jll → Rmath binary" =>
        ("strong copyleft", "GPL-2.0 (R math library); via Makie → StatsFuns → Rmath"),
    "CRlibm_jll → CRlibm binary" => (
        "weak copyleft",
        "CRlibm library LGPL per CRlibm.jl README (binary ships GPL-2.0 text); via Makie",
    ),
    "CRlibm" => ("permissive", "wrapper MIT; wrapped CRlibm library LGPL"),
    "CommonMark" => (
        "permissive",
        "code MIT and BSD-2-Clause; CC-BY-SA-4.0 only for the test spec.txt",
    ),
    "FreeType2_jll → FreeType2 binary" => ("permissive", "dual FTL or GPL-2.0; use FTL"),
    "Zstd_jll → Zstd binary" => ("permissive", "dual BSD-3-Clause or GPL-2.0; use BSD"),
    "XZ_jll → XZ binary" => ("permissive", "liblzma 0BSD; GPL only for helper scripts"),
    "Lz4_jll → Lz4 binary" => ("permissive", "library BSD-2-Clause; GPL only for programs"),
    "MPICH_jll → MPICH binary" => ("permissive", "Argonne permissive license"),
    "Libiconv_jll → Libiconv binary" =>
        ("weak copyleft", "library LGPL, iconv program GPL (GNU libiconv page)"),
    "GettextRuntime_jll → GettextRuntime binary" =>
        ("weak copyleft", "libintl LGPL-2.1; GPL for tools"),
    "Libmount_jll → Libmount binary" =>
        ("weak copyleft", "LGPL-2.1-or-later (util-linux libmount/COPYING)"),
    "Libuuid_jll → Libuuid binary" =>
        ("permissive", "BSD-3-Clause (util-linux libuuid/COPYING)"),
    "Libglvnd_jll → Libglvnd binary" => ("permissive", "MIT-style (README)"),
    "Glib_jll → Glib binary" => ("weak copyleft", "LGPL-2.1-or-later"),
    "FriBidi_jll → FriBidi binary" => ("weak copyleft", "LGPL-2.1-or-later"),
    "LAME_jll → LAME binary" => ("weak copyleft", "LGPL-2.0 (Library GPL)"),
    "Graphite2_jll → Graphite2 binary" =>
        ("permissive", "choice of LGPL-2.1+, MPL-2.0, GPL-2.0+, or MIT; use MIT"),
    "Cairo_jll → Cairo binary" => ("weak copyleft", "LGPL-2.1 or MPL-1.1"),
    "HarfBuzz_jll → HarfBuzz binary" => ("permissive", "Old MIT"),
    "libva_jll → libva binary" => ("permissive", "MIT"),
    "Libtiff_jll → Libtiff binary" => ("permissive", "libtiff license"),
    "Fontconfig_jll → Fontconfig binary" => ("permissive", "HPND-style and MIT"),
    "CoreMath_jll → CoreMath binary" =>
        ("unverified", "no license file in artifact; CoreMath.jl MIT"),
    "Makie → MakieAssets binary" => (
        "review",
        "fonts/LICENSES.md (DejaVu, Noto, Fira, TeX Gyre, Hack); Helvetica.otf has an " *
        "Adobe copyright and no license entry",
    ),
    "HDF5_jll → HDF5 binary" =>
        ("permissive", "BSD-3-Clause (The HDF Group); text not recognized by LicenseCheck"),
    "MicrosoftMPI_jll → MicrosoftMPI binary" => (
        "proprietary, redistributable",
        "Microsoft MPI Redistributable EULA; Windows only; via NCDatasets → NetCDF → HDF5",
    ),
    "aws_c_cal_jll → aws_c_cal binary" => (
        "permissive",
        "Apache-2.0 per upstream awslabs/aws-c-cal LICENSE and NOTICE; via HDF5 ≥ 2",
    ),
    "Git_jll → Git binary" => ("strong copyleft", "GPL-2.0; test only (Documenter)"),
    "OpenSSH_jll → OpenSSH binary" => ("permissive", "BSD-style; test only"),
)
for n in ("X11", "Xau", "Xdmcp", "Xext", "Xfixes", "Xrender")
    REVIEWED["Xorg_lib$(n)_jll → Xorg_lib$(n) binary"] = ("permissive", "HPND/MIT family")
end
REVIEWED["Xorg_xtrans_jll → Xorg_xtrans binary"] = ("permissive", "HPND/MIT family")

function classify(licenses)
    isempty(licenses) && return "unidentified"
    any(l -> occursin(STRONG_COPYLEFT, l), licenses) && return "strong copyleft"
    any(l -> occursin(WEAK_COPYLEFT, l), licenses) && return "weak copyleft"
    all(in(PERMISSIVE), licenses) && return "permissive"
    return "review"
end

# License identifiers found in the license-like files at the top of a package directory.
function scan(dir)
    files, ids = String[], String[]
    isdir(dir) || return files, ids
    for f in find_licenses(dir)
        push!(files, relpath(f.license_filename, dir))
        append!(ids, f.licenses_found)
    end
    return files, unique(ids)
end

# License identifiers found in all files below `dir` (artifact `share/licenses`).
function scan_all(dir)
    files, ids = String[], String[]
    isdir(dir) || return files, ids
    for (root, _, fs) in walkdir(dir), f in fs
        path = joinpath(root, f)
        isfile(path) || continue
        push!(files, relpath(path, dir))
        append!(ids, licensecheck(read(path, String)).licenses_found)
    end
    return files, unique(ids)
end

is_stdlib(info) = info.source !== nothing && startswith(info.source, Sys.STDLIB)

const PLATFORM = let p = get(ENV, "LICENSE_CHECK_PLATFORM", "")
    isempty(p) ? HostPlatform() : parse(Base.BinaryPlatforms.Platform, p)
end

# License files inside the platform artifacts of a `_jll` package.
function artifact_dirs(pkgpath)
    toml = joinpath(pkgpath, "Artifacts.toml")
    isfile(toml) || return Pair{String,String}[]
    out = Pair{String,String}[]
    for name in keys(Pkg.TOML.parsefile(toml))
        meta = artifact_meta(name, toml; platform=PLATFORM)
        meta === nothing && continue
        dir = Pkg.Artifacts.ensure_artifact_installed(name, meta, toml; platform=PLATFORM)
        isdir(dir) && push!(out, name => joinpath(dir, "share", "licenses"))
    end
    return out
end

deps = Pkg.dependencies()
uuid_of = Dict(info.name => uuid for (uuid, info) in deps)
function closure(roots)
    seen = Set{Base.UUID}()
    stack = [uuid_of[r] for r in roots if haskey(uuid_of, r)]
    while !isempty(stack)
        u = pop!(stack)
        u in seen && continue
        push!(seen, u)
        append!(stack, values(deps[u].dependencies))
    end
    return seen
end
runtime = closure(["PlanetTools"])  # Weak dependencies are not part of this closure.
extension = setdiff(closure(WEAKDEPS), runtime)
testonly = setdiff(closure(EXTRAS), runtime, extension)

rows = NamedTuple[]
for (scope, set) in (("runtime", runtime), ("extension", extension), ("test", testonly))
    for u in set
        info = deps[u]
        info.name == "PlanetTools" && continue
        version = info.version === nothing ? "-" : string(info.version)
        if is_stdlib(info) || info.version === nothing
            push!(
                rows,
                (;
                    scope,
                    name=info.name,
                    version,
                    kind="Julia stdlib",
                    files="",
                    ids="MIT; bundled libraries: see Julia THIRDPARTY.md",
                    class="bundled with Julia",
                ),
            )
            continue
        end
        files, ids = scan(info.source)
        push!(
            rows,
            (;
                scope,
                name=info.name,
                version,
                kind=endswith(info.name, "_jll") ? "JLL wrapper" : "Julia package",
                files=join(files, ", "),
                ids=join(ids, ", "),
                class=classify(ids),
            ),
        )
        for (aname, ldir) in artifact_dirs(info.source)
            afiles, aids = scan_all(ldir)
            push!(
                rows,
                (;
                    scope,
                    name="$(info.name) → $aname binary",
                    version,
                    kind="binary artifact",
                    files=join(afiles, ", "),
                    ids=join(aids, ", "),
                    class=isdir(ldir) ? classify(aids) : "no license directory",
                ),
            )
        end
    end
end
sort!(rows; by=r -> (r.scope, lowercase(r.name)))

# Apply the manual review; anything else that is not clearly permissive is flagged.
rows = map(rows) do r
    haskey(REVIEWED, r.name) || return merge(r, (; note=""))
    class, note = REVIEWED[r.name]
    return merge(r, (; class="$class (reviewed)", note))
end
flagged = filter(
    r -> r.class ∉ ("permissive", "bundled with Julia") && isempty(r.note), rows
)
open(OUTFILE, "w") do io
    println(io, "# Third-party licenses of PlanetTools dependencies\n")
    println(
        io,
        "Generated by `scripts/license_check.jl` with Julia $(VERSION) for ",
        "$(triplet(PLATFORM)). Licenses identified with LicenseCheck.jl from the license ",
        "files of each package and of each binary artifact for this platform.\n",
    )
    println(
        io, "| Scope | Package | Version | Kind | Detected license | Class | Review note |"
    )
    println(io, "|:--|:--|:--|:--|:--|:--|:--|")
    for r in rows
        println(
            io,
            "| $(r.scope) | $(r.name) | $(r.version) | $(r.kind) | $(r.ids) | ",
            "$(r.class) | $(replace(r.note, "|" => "\\|")) |",
        )
    end
end

println("Wrote $(length(rows)) entries to $OUTFILE")
for r in flagged
    println(
        "  REVIEW: $(r.name) ($(r.scope)): $(isempty(r.ids) ? "-" : r.ids) [$(r.class)]"
    )
end

# Summary of the classes per scope; also written to the GitHub Actions job summary.
function summary(io)
    println(io, "### License check: ", isempty(flagged) ? "passed" : "needs review", "\n")
    println(io, "| Scope | Class | Entries |\n|:--|:--|--:|")
    for scope in ("runtime", "extension", "test")
        counts = Dict{String,Int}()
        for r in rows
            r.scope == scope && (counts[r.class] = get(counts, r.class, 0) + 1)
        end
        for (class, n) in sort(collect(counts); by=last, rev=true)
            println(io, "| $scope | $class | $n |")
        end
    end
    if !isempty(flagged)
        println(
            io, "\nNot yet reviewed (add to `REVIEWED` in `scripts/license_check.jl`):\n"
        )
        for r in flagged
            println(io, "- $(r.name) ($(r.scope)): $(isempty(r.ids) ? "-" : r.ids)")
        end
    end
    return nothing
end
summary(stdout)
if haskey(ENV, "GITHUB_STEP_SUMMARY")
    open(io -> summary(io), ENV["GITHUB_STEP_SUMMARY"], "a")
end
exit(isempty(flagged) ? 0 : 1)
