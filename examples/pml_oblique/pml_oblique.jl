# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script sweeps the angle of incidence of a plane wave on the truncation of the
periodic vacuum cell of `grazing_80_te.json`, and compares the reflection computed by Palace
with the plane-wave theory.

For each angle θ from 0° to 85°, in TE and TM polarization, the Floquet port at z = 0
launches a plane wave with the Floquet wave vector k_y = k₀ sin θ at 3 GHz. The cell is
terminated by:
1. the first-order absorbing boundary, which reflects (1 - cos θ) / (1 + cos θ) in both
   polarizations;
2. the second-order absorbing boundary, which reflects ((1 - cos θ) / (1 + cos θ))² in TE
   polarization and the same as the first-order one in TM polarization;
3. the meshed PML of `grazing_80_te.json`, for two reflection targets R(0), for which the
   continuous layer reflects R(0)^cos θ;
4. the same PML extruded along the normal of the back of the cell (attribute 8), in the same
   number of elements of the same size, with the periodic boundaries through the layer;
5. the meshed PML on a cell as long as the one in front of the extruded PML, whose mesh
   `mesh/mesh.jl` generates in `postpro/angle/`.

For the absorbing boundaries and the extruded PML, the PML region of the mesh is vacuum and
the truncation is on the back of the cell.

## Prerequisites

This script requires Julia packages. Install them with:

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

## How to run

From the repository root, run:
```bash
julia --project=examples -e 'include("examples/pml_oblique/pml_oblique.jl"); generate_reflection_data()'
```

This requires `palace` to be a runnable command. If it is not, you can pass the path to the
executable, e.g.
```bash
julia --project=examples -e 'include("examples/pml_oblique/pml_oblique.jl"); generate_reflection_data(palace_exec="build/bin/palace")'
```

The script will:
1. Run each case, with its outputs in `postpro/angle`
2. Write the reflection of every case to `postpro/angle/reflection.csv`, and print it with
   the theory and the largest difference between the meshed and the extruded PML
3. Save the figure to `docs/src/assets/examples/pml-reflection-angle.svg`

To redraw the figure from `postpro/angle/reflection.csv` without running Palace, call
`plot_reflection()`.
=#

using CSV
using DataFrames
using JSON
using CairoMakie
using Printf

include(joinpath(@__DIR__, "mesh", "mesh.jl"))

const REPO_ROOT = abspath(joinpath(@__DIR__, "..", ".."))
const FREQUENCY = 3.0e9                       # Hz, as in grazing_80_te.json
const K0 = 2π * FREQUENCY / 299792458.0       # m⁻¹

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
    truncations(targets)

The truncations of the cell, as `(label, R0)` pairs: the absorbing boundaries, labelled
`"Absorbing1"` and `"Absorbing2"` with `R0 = NaN`, and the meshed PML (`"PML"`) and the
extruded PML (`"PMLExtruded"`) for each reflection target `R0`, and the meshed PML on the
cell of the length of that of the extruded PML (`"PMLLong"`).
"""
function truncations(targets)
    return vcat(
        [("Absorbing1", NaN), ("Absorbing2", NaN)],
        [("PML", R0) for R0 in targets],
        [("PMLExtruded", R0) for R0 in targets],
        [("PMLLong", R0) for R0 in targets]
    )
end

"""
    reflection_config(θ, polarization, truncation, R0)

The configuration of one case: the template `grazing_80_te.json` at the angle of incidence
`θ` (degrees), with the `polarization` ("TE" or "TM") and the `truncation` ("Absorbing1",
"Absorbing2", or "PML" or "PMLExtruded" with reflection target `R0`).
"""
function reflection_config(θ, polarization, truncation, R0)
    config = read_config("grazing_80_te.json")
    config["Boundaries"]["Periodic"]["FloquetWaveVector"] = [0.0, K0 * sind(θ), 0.0]
    config["Boundaries"]["FloquetPort"][1]["IncidentPolarization"] = polarization
    if truncation == "PML" || truncation == "PMLLong"
        config["Domains"]["PML"][1]["ReflectionTarget"] = R0
        truncation == "PMLLong" && (config["Model"]["Mesh"] = long_mesh())
    elseif truncation == "PMLExtruded"
        # The layer of the meshed PML: 0.3 m in 24 elements, with the same grading.
        pml = pop!(config["Domains"], "PML")[1]
        delete!(config["Boundaries"], "PEC")
        config["Boundaries"]["PML"] = [
            Dict(
                "Attributes" => [8],
                "Thickness" => pml["Thickness"],
                "Layers" => 24,
                "Order" => pml["Order"],
                "ReflectionTarget" => R0,
                "ReferenceFrequency" => pml["ReferenceFrequency"]
            )
        ]
    else
        delete!(config["Domains"], "PML")
        delete!(config["Boundaries"], "PEC")
        order = truncation == "Absorbing1" ? 1 : 2
        config["Boundaries"]["Absorbing"] = Dict("Attributes" => [8], "Order" => order)
    end
    label =
        isnan(R0) ? lowercase(truncation) : @sprintf("%s_%.0e", lowercase(truncation), R0)
    name = @sprintf("%s_%s_%02d", label, lowercase(polarization), θ)
    config["Problem"]["Output"] = joinpath("postpro", "angle", name)
    return config, name
end

"""
    long_mesh()

The mesh of the cell with 0.5 m in 40 elements between the Floquet port and the meshed PML,
as in front of the PML extruded from the back of the cell of `mesh/pml_oblique.msh`,
generated in `postpro/angle/` if needed.
"""
function long_mesh()
    filename = joinpath(@__DIR__, "postpro", "angle", "pml_oblique_long.msh")
    if !isfile(filename)
        mkpath(dirname(filename))
        generate_pml_oblique_mesh(; filename=filename, L=0.5, n_z=40, verbose=1)
    end
    return filename
end

"""
    run_case(config, name, palace_exec, num_processors)

Run Palace on `config` and return the reflection |S11| in dB of the incident polarization.
"""
function run_case(config, polarization, name, palace_exec, num_processors)
    mkpath(joinpath(@__DIR__, "postpro", "angle"))
    path = joinpath(@__DIR__, "postpro", "angle", "$(name).json")
    open(path, "w") do f
        return JSON.print(f, config, 2)
    end
    run(
        pipeline(
            Cmd(`$palace_exec -np $num_processors $path`; dir=@__DIR__);
            stdout=devnull
        )
    )
    data = CSV.read(
        joinpath(@__DIR__, config["Problem"]["Output"], "port-floquet-S.csv"),
        DataFrame;
        normalizenames=false
    )
    col = findfirst(n -> occursin("|S[P1(0;0)$(polarization)][1]|", n), names(data))
    return data[1, col]
end

"""
    theory_reflection(θ, polarization, truncation, R0)

Plane-wave reflection |R| in dB of the truncation at the angle of incidence `θ` (degrees).
"""
function theory_reflection(θ, polarization, truncation, R0)
    r1 = (1 - cosd(θ)) / (1 + cosd(θ))
    if truncation == "Absorbing1" || (truncation == "Absorbing2" && polarization == "TM")
        return 20 * log10(r1)
    elseif truncation == "Absorbing2"
        return 40 * log10(r1)
    else
        return 20 * cosd(θ) * log10(R0)
    end
end

"""
    generate_reflection_data(;
        palace_exec="palace",
        num_processors::Integer=4,
        angles=0:5:85,
        targets=[1.0e-3, 1.0e-12]
    )

Run every case, write `postpro/angle/reflection.csv`, print the comparison with the theory
and the largest difference between the meshed and the extruded PML, and draw the figure.

# Arguments

  - palace_exec - executable for Palace
  - num_processors - number of processors to use for each simulation
  - angles - angles of incidence, in degrees
  - targets - reflection targets of the PML at normal incidence
"""
function generate_reflection_data(;
    palace_exec="palace",
    num_processors::Integer=4,
    angles=0:5:85,
    targets=[1.0e-3, 1.0e-12]
)
    palace_exec = resolve_executable(palace_exec)
    rows = DataFrame(
        truncation=String[],
        target=Float64[],
        polarization=String[],
        angle_deg=Float64[],
        reflection_dB=Float64[],
        theory_dB=Float64[]
    )
    for (truncation, R0) in truncations(targets), polarization in ("TE", "TM"), θ in angles
        config, name = reflection_config(θ, polarization, truncation, R0)
        R = run_case(config, polarization, name, palace_exec, num_processors)
        Rth = theory_reflection(θ, polarization, truncation, R0)
        push!(rows, (truncation, R0, polarization, θ, R, Rth))
        @printf(
            "%-10s R0 = %-7.0e %s θ = %4.1f°: %9.3f dB (theory %9.3f dB)\n",
            truncation,
            R0,
            polarization,
            θ,
            R,
            Rth
        )
    end
    CSV.write(joinpath(@__DIR__, "postpro", "angle", "reflection.csv"), rows)
    for R0 in targets
        case(t) = rows[(rows.truncation .== t) .& (rows.target .== R0), :reflection_dB]
        extruded = case("PMLExtruded")
        @printf(
            "R0 = %.0e: the extruded PML is within %.3f dB of the meshed PML, and within %.1e dB on a cell of the same length\n",
            R0,
            maximum(abs.(case("PML") .- extruded)),
            maximum(abs.(case("PMLLong") .- extruded))
        )
    end
    plot_reflection()
    return rows
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
const COLORS = Dict(
    "Absorbing1" => RGBf(0.90, 0.62, 0.0),
    "Absorbing2" => RGBf(0.0, 0.62, 0.45),
    "PML1" => RGBf(0.34, 0.71, 0.91),
    "PML2" => RGBf(0.80, 0.47, 0.65)
)

"""
    plot_reflection(;
        csv=joinpath(@__DIR__, "postpro", "angle", "reflection.csv"),
        filename=joinpath(REPO_ROOT, "docs", "src", "assets", "examples",
                          "pml-reflection-angle.svg")
    )

Draw the reflection against the angle of incidence, theory as lines and Palace as markers,
for TE and TM polarization.
"""
function plot_reflection(;
    csv=joinpath(@__DIR__, "postpro", "angle", "reflection.csv"),
    filename=joinpath(
        REPO_ROOT,
        "docs",
        "src",
        "assets",
        "examples",
        "pml-reflection-angle.svg"
    )
)
    data = CSV.read(csv, DataFrame)
    targets = sort(unique(data.target[data.truncation .== "PML"]); rev=true)
    cases = vcat(
        [("Absorbing1", NaN, "First-order absorbing", COLORS["Absorbing1"], :circle)],
        [("Absorbing2", NaN, "Second-order absorbing", COLORS["Absorbing2"], :diamond)],
        [
            (
                "PML",
                R0,
                rich("PML, R(0) = 10", superscript(string(round(Int, log10(R0))))),
                COLORS["PML$(k)"],
                k == 1 ? :rect : :utriangle
            ) for (k, R0) in enumerate(targets)
        ]
    )
    θs = range(0.0, 89.0; length=400)
    fig = with_theme(docs_theme()) do
        fig = Figure(size=(900, 420))
        axes = Axis[]
        for (j, polarization) in enumerate(("TE", "TM"))
            ax = Axis(
                fig[1, j];
                title="$(polarization) polarization",
                xlabel="Angle of incidence (°)",
                ylabel=j == 1 ? "Reflection (dB)" : "",
                xticks=0:15:90,
                limits=(0, 90, -170, 5)
            )
            for (truncation, R0, label, color, marker) in cases
                th = [theory_reflection(θ, polarization, truncation, R0) for θ in θs]
                lines!(ax, θs, max.(th, -1000.0); color=color, linewidth=1.5)
                k =
                    (data.truncation .== truncation) .&
                    (data.polarization .== polarization) .&
                    (isnan(R0) ? trues(nrow(data)) : (data.target .== R0))
                scatter!(
                    ax,
                    data.angle_deg[k],
                    data.reflection_dB[k];
                    color=:transparent,
                    strokecolor=color,
                    strokewidth=1.5,
                    marker=marker,
                    markersize=10
                )
            end
            j == 2 && hideydecorations!(ax; grid=false, ticks=false)
            push!(axes, ax)
        end
        elements = [
            [
                LineElement(; color=color, linewidth=1.5),
                MarkerElement(;
                    color=:transparent,
                    strokecolor=color,
                    strokewidth=1.5,
                    marker=marker,
                    markersize=10
                )
            ] for (_, _, _, color, marker) in cases
        ]
        Legend(
            fig[2, 1:2],
            elements,
            [label for (_, _, label, _, _) in cases];
            orientation=:horizontal,
            nbanks=2,
            tellheight=true,
            framevisible=false
        )
        return fig
    end
    mkpath(dirname(filename))
    save(filename, fig)
    return fig
end
