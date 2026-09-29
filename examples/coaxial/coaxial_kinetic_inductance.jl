# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script runs the magnetostatic coaxial examples and compares the extracted
inductance of a shorted superconducting micro-coax against the exact result.

The script runs:
1. `coaxial_magnetostatic_pec.json`: perfectly conducting inner and outer conductors,
   giving the geometric inductance only
2. `coaxial_magnetostatic_superconductor.json` over a sweep of London penetration depths
   λ, with the film thickness d fixed

For each λ the extracted inductance is compared with

    L = μ₀ℓ/(2π) ln(b/a) + μ₀λ coth(d/λ) ℓ/(2π) (1/a + 1/b)

where the second term is the kinetic inductance of the two superconducting sheets. The
kinetic inductance L - L_PEC is plotted against λ together with the exact sheet
inductance and its thick-film (μ₀λ) and thin-film (μ₀λ²/d) limits.

## Prerequisites

This script requires Julia packages. Install them with:

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

## How to run

From the repository root, run:
```bash
julia --project=examples -e 'include("examples/coaxial/coaxial_kinetic_inductance.jl"); generate_kinetic_inductance_data()'
```

This requires `palace` to be a runnable command. If it is not, you can pass the
path to the executable, e.g.
```bash
julia --project=examples -e 'include("examples/coaxial/coaxial_kinetic_inductance.jl"); generate_kinetic_inductance_data(palace_exec="build/bin/palace")'
```

The script will:
1. Run the PEC case and the superconducting case for each λ
2. Print the extracted and analytic inductances
3. Save a plot of the kinetic inductance against λ to
   `postpro/magnetostatic_kinetic_inductance.png`

## Field images

The field images in the documentation come from the superconducting case on the mesh
refined uniformly twice. With ParaView installed, run from the repository root:
```bash
julia --project=examples -e 'include("examples/coaxial/coaxial_kinetic_inductance.jl"); generate_kinetic_inductance_fields()'
```

This runs the refined case with fields saved to
`postpro/magnetostatic_superconductor_fields`, then calls `pvpython
plot_magnetostatic_docs.py`, which writes `coaxial-5.png` and `coaxial-6.png`. Pass
`pvpython_exec` if `pvpython` is not on your `PATH`.
=#

using CSV
using DataFrames
using JSON
using Measures
using Plots

# Geometry of examples/coaxial/mesh/coaxial.msh with Model.L0 = 1e-6 (μm).
const a_coax = 1.6383e-6 / 2 # Inner conductor radius (m)
const b_coax = 5.461e-6 / 2 # Outer conductor radius (m)
const ℓ_coax = 40.0e-6 # Length (m)
const μ₀ = 4π * 1e-7

"""
    coax_geometric_inductance()

Inductance of a shorted coaxial line with perfectly conducting walls.
"""
coax_geometric_inductance() = μ₀ * ℓ_coax / (2π) * log(b_coax / a_coax)

"""
    coax_kinetic_inductance(λ, d)

Kinetic inductance of the inner and outer conductors, each a superconducting sheet of
thickness `d` and London penetration depth `λ` with field on one side only, so each has
sheet inductance μ₀λ coth(d/λ) per square.
"""
coax_kinetic_inductance(λ, d) =
    μ₀ * λ * coth(d / λ) * ℓ_coax / (2π) * (1 / a_coax + 1 / b_coax)

"""
    extract_inductance(outdir)

Read the self inductance of the single current terminal from `terminal-M.csv`.
"""
function extract_inductance(outdir)
    data = CSV.File(joinpath(outdir, "terminal-M.csv"), header=1) |> DataFrame |> Matrix
    return data[1, 2]
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
    run_config(config, palace_exec, num_processors)

Run Palace on `config` from the example directory, through a temporary configuration file.
"""
function run_config(config, palace_exec, num_processors)
    tmp_path = joinpath(@__DIR__, "coaxial_magnetostatic_tmp.json")
    open(tmp_path, "w") do f
        return JSON.print(f, config, 4)
    end
    run(Cmd(`$palace_exec -np $num_processors $tmp_path`; dir=@__DIR__))
    rm(tmp_path)
    return
end

"""
    generate_kinetic_inductance_data(;
                                      palace_exec::String="palace",
                                      num_processors::Integer=1,
                                      penetration_depths=[0.025, 0.05, 0.1, 0.2, 0.4,
                                                          0.8, 1.6]
                                     )

Generate the data for the superconducting coaxial example and visualize the results.

# Arguments

  - palace_exec - executable for Palace
  - num_processors - number of processors to use for the simulation
  - penetration_depths - London penetration depths to sweep, in μm
"""
function generate_kinetic_inductance_data(;
    palace_exec="palace",
    num_processors::Integer=1,
    penetration_depths=[0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6]
)
    palace_exec = resolve_executable(palace_exec)
    coaxial_dir = @__DIR__

    # Geometric inductance from the PEC case, discarding the terminal output
    run(
        Cmd(
            `$palace_exec -np $num_processors coaxial_magnetostatic_pec.json`;
            dir=coaxial_dir
        )
    )
    L_pec = extract_inductance(joinpath(coaxial_dir, "postpro", "magnetostatic_pec"))

    base_config = read_config("coaxial_magnetostatic_superconductor.json")
    # Thickness in mesh length units, which are μm for this mesh (Model.L0 = 1e-6)
    d = base_config["Boundaries"]["Superconductor"][1]["Thickness"] * 1e-6

    L_sc = similar(penetration_depths, Float64)
    for (i, λ) ∈ enumerate(penetration_depths)
        println("Running superconducting case with λ = $(λ) μm...")
        outdir_rel = "postpro/magnetostatic_superconductor_lambda_$(λ)"
        config = deepcopy(base_config)
        config["Problem"]["Output"] = outdir_rel
        config["Boundaries"]["Superconductor"][1]["PenetrationDepth"] = λ
        run_config(config, palace_exec, num_processors)
        L_sc[i] = extract_inductance(joinpath(coaxial_dir, outdir_rel))
    end

    # Compare with the exact result
    L_geom = coax_geometric_inductance()
    println()
    println(
        "L_PEC = $(round(L_pec * 1e12, digits=4)) pH, " *
        "exact $(round(L_geom * 1e12, digits=4)) pH"
    )
    println("  λ (μm)    L (pH)   exact (pH)   L - L_PEC (pH)   exact kinetic (pH)")
    for (λ, L) ∈ zip(penetration_depths, L_sc)
        L_k = coax_kinetic_inductance(λ * 1e-6, d)
        println(
            lpad(λ, 8),
            lpad(round(L * 1e12, digits=4), 10),
            lpad(round((L_geom + L_k) * 1e12, digits=4), 13),
            lpad(round((L - L_pec) * 1e12, digits=4), 17),
            lpad(round(L_k * 1e12, digits=4), 21)
        )
    end

    # Plot settings
    plotsz = (800, 500)
    fntsz = 12
    fnt = font(fntsz)
    default(
        size=plotsz,
        palette=:Set1_9,
        dpi=300,
        tickfont=fnt,
        guidefont=fnt,
        legendfontsize=fntsz - 2,
        margin=5mm
    )

    λ_fine = exp10.(
        range(
            log10(minimum(penetration_depths) / 2),
            log10(2 * maximum(penetration_depths)),
            length=200
        )
    )
    # Kinetic inductance in pH per unit sheet inductance (H) on the two conductors
    geometric_factor = ℓ_coax / (2π) * (1 / a_coax + 1 / b_coax) * 1e12
    pp = plot(
        xscale=:log10,
        yscale=:log10,
        xlabel="London penetration depth \$\\lambda\$  (μm)",
        ylabel="Kinetic inductance  (pH)",
        legend=:topleft
    )
    plot!(
        pp,
        λ_fine,
        coax_kinetic_inductance.(λ_fine .* 1e-6, d) .* 1e12,
        label="Exact, \$\\mu_0 \\lambda \\coth(d/\\lambda)\$"
    )
    plot!(
        pp,
        λ_fine,
        μ₀ .* λ_fine .* 1e-6 .* geometric_factor,
        label="Thick film, \$\\mu_0 \\lambda\$",
        linestyle=:dash
    )
    plot!(
        pp,
        λ_fine,
        μ₀ .* (λ_fine .* 1e-6) .^ 2 ./ d .* geometric_factor,
        label="Thin film, \$\\mu_0 \\lambda^2 / d\$",
        linestyle=:dash
    )
    scatter!(
        pp,
        penetration_depths,
        (L_sc .- L_pec) .* 1e12,
        label="Palace, \$L - L_{\\mathrm{PEC}}\$",
        color=:black,
        markersize=5
    )

    savefig(pp, joinpath(coaxial_dir, "postpro", "magnetostatic_kinetic_inductance.png"))
    display(pp)

    return
end

"""
    generate_kinetic_inductance_fields(;
                                        palace_exec::String="palace",
                                        num_processors::Integer=1,
                                        pvpython_exec::String="pvpython"
                                       )

Run the superconducting case on the mesh refined uniformly twice, saving the fields, and
render the documentation images with `plot_magnetostatic_docs.py`. The example mesh
resolves the inductance to 0.1% but puts a single element across the gap, so the pointwise
field needs the refinement.

# Arguments

  - palace_exec - executable for Palace
  - num_processors - number of processors to use for the simulation
  - pvpython_exec - executable for ParaView's Python interpreter
"""
function generate_kinetic_inductance_fields(;
    palace_exec="palace",
    num_processors::Integer=1,
    pvpython_exec="pvpython"
)
    palace_exec = resolve_executable(palace_exec)
    pvpython_exec = resolve_executable(pvpython_exec)

    config = read_config("coaxial_magnetostatic_superconductor.json")
    config["Problem"]["Output"] = "postpro/magnetostatic_superconductor_fields"
    config["Model"]["Refinement"] = Dict("UniformLevels" => 2)
    config["Solver"]["Magnetostatic"]["Save"] = 1
    run_config(config, palace_exec, num_processors)

    run(Cmd(`$pvpython_exec plot_magnetostatic_docs.py`; dir=@__DIR__))
    return
end
