# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script runs the patch antenna example and compares its reflection coefficient
S11 with the measured and computed curves of Sheen et al., IEEE Trans. Microw. Theory Tech.
38 (1990) 849-857, Fig. 5, digitized in `sheen1990_s11.csv`.

The script runs `patch_antenna.json`, with the second-order absorbing boundary on the open
faces of the box, and `patch_antenna_pml.json`, with a PML extruded from them, as adaptive
frequency sweeps from 4 to 20 GHz. It prints the frequency and the depth of the resonances of
S11 together with those of the digitized curves, and plots S11 against them.

A resolution study runs the sweep with the absorbing boundary on meshes generated with
`mesh/mesh.jl` for several refinements and finite element orders, and a study of the box
runs both truncations for several distances from the metal to the open faces and sizes of
the wave port.

## Prerequisites

This script requires Julia packages. Install them with:

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

## How to run

From the repository root, run:
```bash
julia --project=examples -e 'include("examples/patch_antenna/patch_antenna.jl"); generate_patch_antenna_data()'
```

This requires `palace` to be a runnable command. If it is not, you can pass the path to the
executable, e.g.
```bash
julia --project=examples -e 'include("examples/patch_antenna/patch_antenna.jl"); generate_patch_antenna_data(palace_exec="build/bin/palace")'
```

The script will:
1. Run each case, with its outputs in `postpro/`
2. Print the frequency and the depth of the resonances
3. Save the figure of S11 to `docs/src/assets/examples/patch_antenna-s11.svg`

The resolution study runs with
```bash
julia --project=examples -e 'include("examples/patch_antenna/patch_antenna.jl"); generate_resolution_study()'
```
and writes its meshes and outputs to `postpro/resolution/`, and the study of the box with
```bash
julia --project=examples -e 'include("examples/patch_antenna/patch_antenna.jl"); generate_box_study()'
```
to `postpro/box/`. After `generate_patch_antenna_data()`, the field plot of the first
resonance runs with
```bash
julia --project=examples -e 'include("examples/patch_antenna/patch_antenna.jl"); generate_field_plot()'
```
and saves `docs/src/assets/examples/patch_antenna-field.png`.
=#

using CSV
using DataFrames
using JSON
using CairoMakie
using Printf

include(joinpath(@__DIR__, "mesh", "mesh.jl"))

const REPO_ROOT = abspath(joinpath(@__DIR__, "..", ".."))

# Frequency windows (GHz) of the resonances of S11 below 20 GHz.
const RESONANCE_WINDOWS =
    [(5.5, 6.8), (7.0, 8.0), (9.3, 10.6), (11.4, 12.8), (13.8, 15.0), (17.0, 18.9)]

"""
    read_config(file)

Read a configuration file from the example directory, stripping the // comments that
Palace accepts but JSON.jl does not.
"""
function read_config(file)
    return open(joinpath(@__DIR__, file)) do f
        return JSON.parse(replace(read(f, String), r"//[^\n]*" => ""))
    end
end

"""
    resolve_executable(exec)

Return `exec` as an absolute path if it is a relative path, since the simulations run from
the example directory. Bare command names are left for `PATH` lookup.
"""
function resolve_executable(exec)
    is_path = occursin(Base.Filesystem.path_separator, exec)
    return is_path && !isabspath(exec) ? abspath(exec) : exec
end

"""
    run_config(config, palace_exec, num_processors)

Run Palace on `config` from the example directory, through a configuration file written
next to its outputs.
"""
function run_config(config, palace_exec, num_processors)
    outdir = joinpath(@__DIR__, config["Problem"]["Output"])
    mkpath(outdir)
    path = outdir * ".json"
    open(path, "w") do f
        return JSON.print(f, config, 2)
    end
    run(Cmd(`$palace_exec -np $num_processors $path`; dir=@__DIR__))
    return outdir
end

"""
    read_s11(outdir)

Read the frequency (GHz) and |S11| (dB) from `port-S.csv` in `outdir`.
"""
function read_s11(outdir)
    data = CSV.read(joinpath(outdir, "port-S.csv"), DataFrame; normalizenames=false)
    return data[!, 1], data[!, 2]
end

"""
    read_sheen()

Read the digitized measured and computed curves of Sheen et al. (1990), Fig. 5, as
`(f_meas, s_meas, f_calc, s_calc)`.
"""
function read_sheen()
    data = CSV.read(joinpath(@__DIR__, "sheen1990_s11.csv"), DataFrame; comment="#")
    meas = data[data.curve .== "measurement", :]
    calc = data[data.curve .== "calculation", :]
    return meas.frequency_GHz, meas.s11_dB, calc.frequency_GHz, calc.s11_dB
end

"""
    resonances(f, s)

Frequency and depth of the minimum of |S11| `s` (dB) in each of the `RESONANCE_WINDOWS`,
refined by a parabola through the sample at the minimum and its neighbors.
"""
function resonances(f, s)
    out = Tuple{Float64, Float64}[]
    for (lo, hi) in RESONANCE_WINDOWS
        k = findall(x -> lo <= x <= hi, f)
        i = k[argmin(s[k])]
        fmin, smin = f[i], s[i]
        if 1 < i < length(f) && abs(f[i + 1] - f[i] - (f[i] - f[i - 1])) < 1.0e-9
            y0, y1, y2 = s[i - 1], s[i], s[i + 1]
            den = y0 - 2y1 + y2
            if den > 0
                δ = (y0 - y2) / 2den
                fmin += δ * (f[i + 1] - f[i])
                smin = y1 - (y0 - y2) * δ / 4
            end
        end
        push!(out, (fmin, smin))
    end
    return out
end

"""
    print_resonances(cases)

Print the resonances of each `(label, f, s)` in `cases` next to those of the digitized curves.
"""
function print_resonances(cases)
    fm, sm, fc, sc = read_sheen()
    all = vcat([("Measured", fm, sm), ("Sheen et al., FDTD", fc, sc)], cases)
    println("\nResonances of S11 (GHz / dB):")
    for (label, f, s) in all
        @printf("  %-42s", label)
        for (fr, sr) in resonances(f, s)
            @printf(" %7.3f / %6.2f", fr, sr)
        end
        println()
    end
    return
end

"""
    truncation_cases()

The truncations of the open faces compared by the script, as `(name, label, file)`: the
second-order absorbing boundary of `patch_antenna.json`, and the extruded PML of
`patch_antenna_pml.json`.
"""
function truncation_cases()
    return [
        ("abc2", "second-order absorbing boundary", "patch_antenna.json"),
        ("pml", "PML", "patch_antenna_pml.json")
    ]
end

"""
    generate_patch_antenna_data(;
        palace_exec="palace",
        num_processors::Integer=16,
        order::Integer=2
    )

Run the adaptive sweep of `patch_antenna.json` for each of the `truncation_cases()`, print
their resonances, and plot S11.

# Arguments

  - palace_exec - executable for Palace
  - num_processors - number of processors to use for the simulation
  - order - finite element order
"""
function generate_patch_antenna_data(;
    palace_exec="palace",
    num_processors::Integer=16,
    order::Integer=2
)
    palace_exec = resolve_executable(palace_exec)
    results, plotted = [], []
    for (name, label, file) in truncation_cases()
        config = read_config(file)
        config["Solver"]["Order"] = order
        config["Problem"]["Output"] = joinpath("postpro", name)
        t = @elapsed outdir = run_config(config, palace_exec, num_processors)
        @printf("%s: %.0f s\n", label, t)
        push!(results, ("Palace, " * label, read_s11(outdir)...))
        push!(plotted, (label, outdir))
    end
    print_resonances(results)
    plot_s11(plotted)
    return
end

"""
    run_mesh_study(name, cases, palace_exec, num_processors; truncations=truncation_cases())

Run the adaptive sweep of each of the `truncations` for each `(mesh_kwargs, order)` in
`cases`, on a mesh generated by `generate_patch_antenna_mesh(; mesh_kwargs...)`, with
outputs in `postpro/name/`, and print the resonances of S11 of each case.
"""
function run_mesh_study(
    name,
    cases,
    palace_exec,
    num_processors;
    truncations=truncation_cases()
)
    palace_exec = resolve_executable(palace_exec)
    dir = joinpath(@__DIR__, "postpro", name)
    mkpath(dir)
    results = []
    for (mesh_kwargs, order) in cases, (trunc, _, file) in truncations
        tag = join(["$(k)_$(v)" for (k, v) in pairs(mesh_kwargs)], "_")
        mesh = joinpath(dir, "$(tag).msh")
        isfile(mesh) ||
            generate_patch_antenna_mesh(; mesh_kwargs..., filename=mesh, verbose=1)
        config = read_config(file)
        config["Model"]["Mesh"] = mesh
        config["Solver"]["Order"] = order
        config["Problem"]["Output"] =
            joinpath("postpro", name, "$(tag)_order_$(order)_$(trunc)")
        t = @elapsed outdir = run_config(config, palace_exec, num_processors)
        f, s = read_s11(outdir)
        dofs =
            JSON.parsefile(joinpath(outdir, "palace.json"))["Problem"]["DegreesOfFreedom"]
        push!(
            results,
            (
                @sprintf(
                    "%s, order %d, %s (%d unknowns, %.0f s)",
                    tag,
                    order,
                    trunc,
                    dofs,
                    t
                ),
                f,
                s
            )
        )
    end
    print_resonances(results)
    return
end

"""
    generate_resolution_study(;
        palace_exec="palace",
        num_processors::Integer=16,
        refinements=[1.0, 1.4, 2.0],
        orders=[2, 3]
    )

Run the adaptive sweep of `patch_antenna.json` on meshes generated with `mesh/mesh.jl` for
each of the `refinements` (which divide the element sizes), with each of the finite element
`orders`, and print the resonances of S11 of each case.
"""
function generate_resolution_study(;
    palace_exec="palace",
    num_processors::Integer=16,
    refinements=[1.0, 1.4, 2.0],
    orders=[2, 3]
)
    cases = [((refinement=r,), p) for r in refinements for p in orders]
    return run_mesh_study(
        "resolution",
        cases,
        palace_exec,
        num_processors;
        truncations=truncation_cases()[1:1]
    )
end

"""
    generate_box_study(;
        palace_exec="palace",
        num_processors::Integer=16,
        margins=[8.0, 10.0, 15.0],
        ports=[(12.0, 4.0)],
        edges=[4.0],
        order::Integer=2
    )

Run the adaptive sweep of each of the `truncation_cases()` on meshes generated with
`mesh/mesh.jl` with other distances from the metal to the open faces (`margins`, in mm),
other sizes of the wave port (`ports`, as `(width, height)` in mm), and other element sizes
at the edges between the open faces (`edges`, in mm), and print the resonances of S11 of
each case.
"""
function generate_box_study(;
    palace_exec="palace",
    num_processors::Integer=16,
    margins=[8.0, 10.0, 15.0],
    ports=[(12.0, 4.0)],
    edges=[4.0],
    order::Integer=2
)
    cases = vcat(
        [((margin=m,), order) for m in margins],
        [((port_width=w, port_height=h), order) for (w, h) in ports],
        [((size_open_edge=e,), order) for e in edges]
    )
    return run_mesh_study("box", cases, palace_exec, num_processors)
end

"""
    generate_field_plot(;
        palace_exec="palace",
        num_processors::Integer=16,
        order::Integer=2,
        frequency=nothing,
        spacing=0.2,
        filename=joinpath(REPO_ROOT, "docs", "src", "assets", "examples",
                          "patch_antenna-field.png")
    )

Plot the normal electric field in the substrate, a quarter of its thickness above the ground plane, at the first resonance of
the patch, with the PML of `patch_antenna_pml.json` (by default at the frequency of the
minimum of S11 near 7.5 GHz in `postpro/pml`), sampled with a grid of probes of the given `spacing` (mm).
"""
function generate_field_plot(;
    palace_exec="palace",
    num_processors::Integer=16,
    order::Integer=2,
    frequency=nothing,
    spacing=0.2,
    filename=joinpath(
        REPO_ROOT,
        "docs",
        "src",
        "assets",
        "examples",
        "patch_antenna-field.png"
    )
)
    palace_exec = resolve_executable(palace_exec)
    if isnothing(frequency)
        f, s = read_s11(joinpath(@__DIR__, "postpro", "pml"))
        frequency = resonances(f, s)[2][1]
    end
    h, L1, L2, y_patch = 0.794, 12.45, 16.0, 20.0
    xs = collect((-L1 / 2 - 2.0):spacing:(L1 / 2 + 2.0))
    ys = collect((y_patch - 4.0):spacing:(y_patch + L2 + 2.0))
    config = read_config("patch_antenna_pml.json")
    config["Solver"]["Order"] = order
    config["Solver"]["Driven"] =
        Dict("Samples" => [Dict("Type" => "Point", "Freq" => [frequency])], "SaveStep" => 0)
    config["Problem"]["Output"] = joinpath("postpro", "field")
    probes = [(x, y) for y in ys for x in xs]
    config["Domains"]["Postprocessing"] = Dict(
        "Probe" => [
            Dict("Index" => i, "Center" => [x, y, h / 4]) for
            (i, (x, y)) in enumerate(probes)
        ]
    )
    outdir = run_config(config, palace_exec, num_processors)
    data = CSV.read(joinpath(outdir, "probe-E.csv"), DataFrame; normalizenames=false)
    index = Dict(strip(n) => k for (k, n) in enumerate(names(data)))
    col(part, i) = data[1, index["$(part){E_z[$(i)]} (V/m)"]]
    ez = [complex(col("Re", i), col("Im", i)) for i in eachindex(probes)]
    # Rotate the phase so that the field is real where it is largest.
    k = argmax(abs.(ez))
    ez .*= conj(ez[k]) / abs(ez[k])
    Ez = reshape(real.(ez), length(xs), length(ys)) ./ abs(ez[k])
    fig = with_theme(docs_theme()) do
        fig = Figure(size=(560, 640))
        ax = Axis(
            fig[1, 1];
            xlabel="x (mm)",
            ylabel="y (mm)",
            aspect=DataAspect(),
            title=rich(
                "Re{E",
                subscript("z"),
                @sprintf("} at %.3f GHz, z = h / 4", frequency)
            )
        )
        hm = heatmap!(ax, xs, ys, Ez; colormap=:RdBu, colorrange=(-1, 1), interpolate=true)
        x_line = -L1 / 2 + 2.09
        lines!(
            ax,
            [
                -L1 / 2,
                L1 / 2,
                L1 / 2,
                -L1 / 2,
                -L1 / 2,
                x_line,
                x_line,
                NaN,
                x_line + 2.46,
                x_line + 2.46
            ],
            [
                y_patch,
                y_patch,
                y_patch + L2,
                y_patch + L2,
                y_patch,
                y_patch,
                ys[1],
                NaN,
                ys[1],
                y_patch
            ];
            color=RGBf(0.5, 0.5, 0.5),
            linewidth=1.5
        )
        Colorbar(fig[1, 2], hm; label="Normalized")
        return fig
    end
    mkpath(dirname(filename))
    save(filename, fig; px_per_unit=2)
    return frequency
end

"""
    docs_theme()

A theme for figures which read on the light and the dark theme of the documentation: a
transparent background with gray text and axes.
"""
function docs_theme()
    gray = RGBf(0.5, 0.5, 0.5)
    return Theme(
        backgroundcolor=:transparent,
        textcolor=gray,
        fontsize=16,
        Axis=(
            backgroundcolor=:transparent,
            bottomspinecolor=gray,
            leftspinecolor=gray,
            topspinecolor=gray,
            rightspinecolor=gray,
            xtickcolor=gray,
            ytickcolor=gray,
            xgridcolor=RGBAf(0.5, 0.5, 0.5, 0.25),
            ygridcolor=RGBAf(0.5, 0.5, 0.5, 0.25)
        ),
        Legend=(backgroundcolor=:transparent, framecolor=gray)
    )
end

# Colors which read on light and dark backgrounds (Okabe-Ito).
const COLORS = [RGBf(0.0, 0.62, 0.45), RGBf(0.34, 0.71, 0.91), RGBf(0.80, 0.47, 0.65)]

"""
    plot_s11(cases;
        filename=joinpath(REPO_ROOT, "docs", "src", "assets", "examples",
                          "patch_antenna-s11.svg")
    )

Plot |S11| of each `(label, outdir)` in `cases` on the digitized measured and computed
curves of Sheen et al. (1990): from 4 to 20 GHz, and around the first resonance.
"""
function plot_s11(
    cases;
    filename=joinpath(
        REPO_ROOT,
        "docs",
        "src",
        "assets",
        "examples",
        "patch_antenna-s11.svg"
    )
)
    fm, sm, fc, sc = read_sheen()
    gray = RGBf(0.5, 0.5, 0.5)
    fig = with_theme(docs_theme()) do
        fig = Figure(size=(900, 640))
        ax1 = Axis(
            fig[1, 1:2];
            xlabel="Frequency (GHz)",
            ylabel="|S₁₁| (dB)",
            xticks=4:2:20,
            limits=(4, 20, -40, 1)
        )
        ax2 = Axis(
            fig[2, 1];
            xlabel="Frequency (GHz)",
            ylabel="|S₁₁| (dB)",
            limits=(7.0, 8.0, -30, 1)
        )
        ax3 = Axis(fig[2, 2]; xlabel="Frequency (GHz)", limits=(17.0, 19.0, -40, 1))
        for ax in (ax1, ax2, ax3)
            lines!(ax, fm, sm; color=gray, linestyle=:dash, linewidth=1.5, label="Measured")
            lines!(ax, fc, sc; color=gray, linewidth=1.5, label="Sheen et al., FDTD")
            for (k, (label, outdir)) in enumerate(cases)
                f, s = read_s11(outdir)
                lines!(ax, f, s; color=COLORS[k], linewidth=2, label="Palace, " * label)
            end
        end
        Legend(fig[3, 1:2], ax1; orientation=:horizontal, framevisible=false)
        return fig
    end
    mkpath(dirname(filename))
    save(filename, fig)
    return fig
end
