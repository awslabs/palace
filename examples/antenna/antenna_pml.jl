# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script compares the truncations of the outer sphere of the half-wave dipole
antenna: the first- and second-order absorbing boundary of `antenna_halfwave_dipole.json`,
and the PML extruded from the sphere of `antenna_halfwave_dipole_pml.json` for several
numbers of layers and thicknesses.

For each case, it prints the input impedance of the antenna, Z = R (1 + S11) / (1 - S11)
for the port resistance R = 50 Ω, and the directivity of the far field computed on the
outer sphere (attribute 4, an interior boundary when the PML is extruded from it). It then
runs the second-order absorbing boundary and the PML of `antenna_halfwave_dipole_pml.json`
on a finer mesh, and plots their E-plane pattern against the pattern of the thin half-wave
dipole, |cos(π/2 cos θ) / sin θ|, whose directivity is 1.64.

## Prerequisites

This script requires Julia packages. Install them with:

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

## How to run

From the repository root, run:
```bash
julia --project=examples -e 'include("examples/antenna/antenna_pml.jl"); generate_antenna_pml_data()'
```

This requires `palace` to be a runnable command. If it is not, you can pass the path to the
executable, e.g.
```bash
julia --project=examples -e 'include("examples/antenna/antenna_pml.jl"); generate_antenna_pml_data(palace_exec="build/bin/palace")'
```

The script will:
1. Run each case, with its outputs in `postpro/truncation/`, where it also generates the
   finer mesh
2. Print the input impedance and the directivity of each case
3. Save the figure of the E-plane pattern to
   `docs/src/assets/examples/antenna-pml-eplane.svg`

To redraw the figure from the outputs without running Palace, call `plot_eplane()`.
=#

include(joinpath(@__DIR__, "plot_farfield.jl"))
include(joinpath(@__DIR__, "mesh", "mesh.jl"))

import JSON
import Printf: @sprintf

const REPO_ROOT = abspath(joinpath(@__DIR__, "..", ".."))

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
    truncation_cases(; thicknesses=[0.5, 1.0, 2.0], layers=[2, 4, 8])

The truncations compared by `generate_antenna_pml_data`, as `(name, label, config)`: the
first- and second-order absorbing boundary, and the extruded PML with each of the
`thicknesses` (m) and numbers of `layers`.
"""
function truncation_cases(; thicknesses=[0.5, 1.0, 2.0], layers=[2, 4, 8])
    cases = []
    for order in (1, 2)
        config = read_config("antenna_halfwave_dipole.json")
        config["Boundaries"]["Absorbing"]["Order"] = order
        label = order == 1 ? "First-order absorbing" : "Second-order absorbing"
        push!(cases, ("absorbing_order_$(order)", label, config))
    end
    for t in thicknesses, n in layers
        config = read_config("antenna_halfwave_dipole_pml.json")
        config["Boundaries"]["PML"][1]["Thickness"] = t
        config["Boundaries"]["PML"][1]["Layers"] = n
        push!(cases, ("pml_$(t)_$(n)", @sprintf("PML, %.1f m, %d layers", t, n), config))
    end
    return cases
end

"""
    read_impedance(outdir; R=50.0)

The input impedance `R (1 + S11) / (1 - S11)` of the port, from `port-S.csv` in `outdir`.
"""
function read_impedance(outdir; R=50.0)
    data = CSV.read(joinpath(outdir, "port-S.csv"), DataFrame; normalizenames=false)
    S = 10^(data[1, 2] / 20) * cis(deg2rad(data[1, 3]))
    return R * (1 + S) / (1 - S)
end

"""
    run_case(name, label, config, palace_exec, num_processors, nsample; mesh=nothing)

Run `config` with its outputs in `postpro/truncation/name`, `nsample` points on the far-field
sphere, the complex coarse solve and, if given, the `mesh`, and return `(label, Z, D, t)`: the
input impedance, the directivity and the wall time.
"""
function run_case(name, label, config, palace_exec, num_processors, nsample; mesh=nothing)
    outdir = joinpath("postpro", "truncation", name)
    config["Problem"]["Output"] = outdir
    config["Problem"]["OutputFormats"] = Dict("Paraview" => false)
    config["Boundaries"]["Postprocessing"] =
        Dict("FarField" => Dict("Attributes" => [4], "NSample" => nsample))
    # The domain is three wavelengths across, where the complex coarse solve needs far fewer
    # iterations (the real one does not converge on the finer mesh).
    config["Solver"]["Linear"]["ComplexCoarseSolve"] = true
    isnothing(mesh) || (config["Model"]["Mesh"] = mesh)
    path = joinpath(@__DIR__, outdir * ".json")
    mkpath(dirname(path))
    open(path, "w") do f
        return JSON.print(f, config, 2)
    end
    log = joinpath(@__DIR__, outdir * ".log")
    t = @elapsed run(
        pipeline(Cmd(`$palace_exec -np $num_processors $path`; dir=@__DIR__); stdout=log)
    )
    occursin("did NOT converge", read(log, String)) &&
        error("The linear solver did not converge for $(name), see $(log)")
    outdir = joinpath(@__DIR__, outdir)
    data, _ = read_farfield(joinpath(outdir, "farfield-rE.csv"))
    return (label, read_impedance(outdir), directivity(data), t)
end

"""
    print_results(results)

Print the input impedance and the directivity of each `(label, Z, D, t)` in `results`, and
the distance of the impedance from that of the last case.
"""
function print_results(results)
    Z_ref = results[end][2]
    for (label, Z, D, t) in results
        @printf(
            "  %-34s Z = %7.2f %s %6.2fj Ω, |Z - Z_ref| = %5.2f Ω, D = %.3f (%4.0f s)\n",
            label,
            real(Z),
            imag(Z) < 0 ? "-" : "+",
            abs(imag(Z)),
            abs(Z - Z_ref),
            D,
            t
        )
    end
    println("  (Z_ref is the last case)")
    return
end

"""
    generate_antenna_pml_data(;
        palace_exec="palace",
        num_processors::Integer=16,
        nsample::Integer=64800,
        thicknesses=[0.5, 1.0, 2.0],
        layers=[2, 4, 8],
        n_farfield=10
    )

Run each of the `truncation_cases` on `mesh/antenna.msh`, with `nsample` points on the
far-field sphere, and print the input impedance and the directivity of each. Then run the
second-order absorbing boundary and the PML of `antenna_halfwave_dipole_pml.json` on a finer
mesh, with `n_farfield` elements per wavelength away from the antenna (generated in
`postpro/truncation/`), print the same, and plot their E-plane pattern.

# Arguments

  - palace_exec - executable for Palace
  - num_processors - number of processors to use for the simulation
  - nsample - number of points on the far-field sphere
  - thicknesses - thicknesses of the PML (m)
  - layers - numbers of layers of the PML
  - n_farfield - number of elements per wavelength away from the antenna of the finer mesh
"""
function generate_antenna_pml_data(;
    palace_exec="palace",
    num_processors::Integer=16,
    nsample::Integer=64800,
    thicknesses=[0.5, 1.0, 2.0],
    layers=[2, 4, 8],
    n_farfield=10
)
    palace_exec = resolve_executable(palace_exec)
    results = [
        run_case(name, label, config, palace_exec, num_processors, nsample) for
        (name, label, config) in truncation_cases(; thicknesses=thicknesses, layers=layers)
    ]
    println("\nInput impedance and directivity on mesh/antenna.msh:")
    print_results(results)

    mesh = joinpath(@__DIR__, "postpro", "truncation", "antenna_$(n_farfield).msh")
    isfile(mesh) || generate_antenna_mesh(; filename=mesh, n_farfield=n_farfield, verbose=1)
    fine = [
        run_case(
            "fine_" * name,
            label * ", finer mesh",
            config,
            palace_exec,
            num_processors,
            nsample;
            mesh=mesh
        ) for
        (name, label, config) in truncation_cases(; thicknesses=[1.0], layers=[4]) if
        name != "absorbing_order_1"
    ]
    println("\nInput impedance and directivity with $(n_farfield) elements per wavelength:")
    print_results(fine)
    data, _ = read_farfield(
        joinpath(@__DIR__, "postpro", "truncation", "fine_pml_1.0_4", "farfield-rE.csv")
    )
    @printf(
        "Thin half-wave dipole on the same points: D = %.3f\n",
        directivity(data; magnitude=halfwave_dipole_pattern.(data[:, "theta"]))
    )
    plot_eplane()
    return results, fine
end

"""
    docs_theme()

A theme for figures which read on the light and the dark theme of the documentation: a
transparent background with gray text and axes.
"""
function docs_theme()
    gray = Makie.RGBf(0.5, 0.5, 0.5)
    return Makie.Theme(
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
            xgridcolor=Makie.RGBAf(0.5, 0.5, 0.5, 0.25),
            ygridcolor=Makie.RGBAf(0.5, 0.5, 0.5, 0.25)
        ),
        Legend=(backgroundcolor=:transparent, framecolor=gray)
    )
end

"""
    eplane_cut(data)

The far field `|r E|` on the E-plane (the xz-plane: φ = 0° and 180°), normalized to its
maximum, as `(θ, E)` with θ from 0° to 180° on both halves of the plane.
"""
function eplane_cut(data)
    mask = (abs.(data[:, "phi"]) .< 1.0e-6) .| (abs.(data[:, "phi"] .- 180) .< 1.0e-6)
    E = compute_field_magnitude(data)
    θ, Ecut = data[mask, "theta"], E[mask]
    p = sortperm(θ)
    return θ[p], Ecut[p] ./ maximum(E)
end

"""
    plot_eplane(;
        cases=[
            ("fine_absorbing_order_2", "Second-order absorbing"),
            ("fine_pml_1.0_4", "PML, 1.0 m, 4 layers")
        ],
        filename=joinpath(REPO_ROOT, "docs", "src", "assets", "examples",
                          "antenna-pml-eplane.svg")
    )

Plot the E-plane pattern of each `(name, label)` in `cases`, from the outputs in
`postpro/truncation/name`, against the pattern of the thin half-wave dipole, and below it
the difference from the thin dipole in dB. Print the largest difference from the thin dipole,
relative to the maximum, and the largest difference between the cases.
"""
function plot_eplane(;
    cases=[
        ("fine_absorbing_order_2", "Second-order absorbing"),
        ("fine_pml_1.0_4", "PML, 1.0 m, 4 layers")
    ],
    filename=joinpath(
        REPO_ROOT,
        "docs",
        "src",
        "assets",
        "examples",
        "antenna-pml-eplane.svg"
    )
)
    gray = Makie.RGBf(0.5, 0.5, 0.5)
    colors = [Makie.RGBf(0.34, 0.71, 0.91), Makie.RGBf(0.0, 0.62, 0.45)]
    θth = range(0, 180, length=721)
    Eth = halfwave_dipole_pattern.(θth)
    cuts = [
        eplane_cut(
            first(
                read_farfield(
                    joinpath(@__DIR__, "postpro", "truncation", name, "farfield-rE.csv")
                )
            )
        ) for (name, _) in cases
    ]
    for ((_, label), (θ, E)) in zip(cases, cuts)
        @printf(
            "%s: largest difference from the thin dipole %.4f of the maximum\n",
            label,
            maximum(abs.(E .- halfwave_dipole_pattern.(θ)))
        )
    end
    @printf(
        "Largest difference between the cases: %.2e of the maximum\n",
        maximum(abs.(cuts[1][2] .- cuts[end][2]))
    )
    fig = Makie.with_theme(docs_theme()) do
        fig = Makie.Figure(size=(800, 640))
        ax1 = Makie.Axis(
            fig[1, 1];
            ylabel="|E| / max |E| (dB)",
            xticks=0:30:180,
            limits=(0, 180, -30, 1)
        )
        ax2 = Makie.Axis(
            fig[2, 1];
            xlabel="θ (deg.)",
            ylabel="Difference from the thin dipole (dB)",
            xticks=0:30:180,
            limits=(0, 180, nothing, nothing)
        )
        Makie.hidexdecorations!(ax1; grid=false, ticks=false)
        Makie.lines!(
            ax1,
            θth,
            20 .* log10.(max.(Eth, 1e-3));
            color=gray,
            linewidth=1.5,
            label="Thin half-wave dipole"
        )
        for (k, ((_, label), (θ, E))) in enumerate(zip(cases, cuts))
            Makie.scatter!(
                ax1,
                θ[1:4:end],
                20 .* log10.(E[1:4:end]);
                color=colors[k],
                markersize=k == 1 ? 9 : 5,
                label="Palace, " * label
            )
            interior = 1.0 .< θ .< 179.0
            Makie.lines!(
                ax2,
                θ[interior],
                20 .* log10.(E[interior] ./ halfwave_dipole_pattern.(θ[interior]));
                color=colors[k],
                linewidth=k == 1 ? 4 : 1.5
            )
        end
        Makie.rowsize!(fig.layout, 2, Makie.Relative(0.35))
        Makie.Legend(fig[3, 1], ax1; orientation=:horizontal, framevisible=false)
        return fig
    end
    mkpath(dirname(filename))
    Makie.save(filename, fig)
    return fig
end
