# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script runs and verifies the rational-impedance boundary example.

It models a short section of a parallel-plate (TEM) line terminated, at port 2, by a
parallel-RLC sheet. The same termination is realized two ways and the reflection S11 at the
source port (port 1) is compared:

  1. `parallel_rlc_lumped.json`   - port 2 is a `LumpedPort` with surface Rs/Ls/Cs (reference)
  2. `parallel_rlc_rational.json` - port 2 is the new `RationalImpedance` boundary

Geometry (mesh/mesh.jl): plates y=0,y=h are PEC, side walls x=0,x=w are PMC (clean TEM mode),
z=0 is the excited source port, z=L is the termination under test. With Cartesian uniform
lumped ports the per-square admittance applied by the `LumpedPort` equals exactly that of the
`RationalImpedance`, so the two S11 curves must coincide to solver tolerance.

The termination is a parallel RLC tank resonant at 10 GHz with Rs = η (matched at resonance),
giving a clear, strongly frequency-dependent |S11| notch at 10 GHz. The rational coefficients
(highest-degree-first, s = iω in rad/s) follow from Zs = 1/(1/Rs + 1/(sLs) + sCs):

  Numerator   = [Ls·Rs, 0]
  Denominator = [Ls·Rs·Cs, Ls, Rs]

## Prerequisites

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

## How to run

Generate the mesh (once):
```bash
julia --project=examples -e 'include("examples/rational_impedance/mesh/mesh.jl"); generate_parallel_plate_mesh(filename="parallel_plate.msh")'
```

Run both simulations and compare:
```bash
julia --project=examples -e 'include("examples/rational_impedance/rational_impedance.jl"); generate_rational_impedance_data(num_processors=4)'
```

If `palace` is not on PATH, pass the executable:
```bash
... generate_rational_impedance_data(palace_exec="build/bin/palace", num_processors=4)'
```

## Output

Prints the maximum difference between the two S11 curves (PASS if below tolerance) and saves
an overlay plot to `postpro/rational_impedance.png`.
=#

using CSV
using DataFrames
using Measures
using Plots

"""
    generate_rational_impedance_data(; palace_exec="palace", num_processors=1)

Run the `LumpedPort` (reference) and `RationalImpedance` (under test) configurations, compare
the reflection S11 at the source port, and plot the result.

# Arguments

  - palace_exec - executable for Palace
  - num_processors - number of processors to use for the simulation
"""
function generate_rational_impedance_data(; palace_exec="palace", num_processors::Integer=1)
    palace_exec_is_path = occursin(Base.Filesystem.path_separator, palace_exec)
    if palace_exec_is_path
        palace_exec = isabspath(palace_exec) ? palace_exec : abspath(palace_exec)
    end

    example_dir = @__DIR__

    # Run both configurations.
    for cfg ∈ ["parallel_rlc_lumped.json", "parallel_rlc_rational.json"]
        base_cmd = `$palace_exec -np $num_processors $cfg`
        run(Cmd(base_cmd; dir=example_dir))
    end

    # Parse S-parameters. Select the |S[1][1]| (dB) and arg(S[1][1]) (deg.) columns by name,
    # since the reference run has extra (port-2) columns.
    function load_s11(subdir)
        file = joinpath(example_dir, "postpro", subdir, "port-S.csv")
        df = CSV.File(file, header=1) |> DataFrame
        names_ = names(df)
        fcol = names_[findfirst(n -> startswith(strip(n), "f "), names_)]
        mcol = names_[findfirst(n -> occursin("|S[1][1]|", n), names_)]
        acol = names_[findfirst(n -> occursin("arg(S[1][1])", n), names_)]
        return Float64.(df[!, fcol]), Float64.(df[!, mcol]), Float64.(df[!, acol])
    end

    f, mag_ref, ang_ref = load_s11("lumped")
    _, mag_rat, ang_rat = load_s11("rational")

    dmag = maximum(abs.(mag_ref .- mag_rat))
    dang = maximum(abs.(ang_ref .- ang_rat))
    println("\nmax |Δ|S11||     = $(dmag) dB")
    println("max |Δarg(S11)|  = $(dang) deg")
    println(
        (dmag < 1e-4 && dang < 1e-2) ? "PASS: RationalImpedance matches LumpedPort RLC" :
        "CHECK: discrepancy exceeds tolerance"
    )

    # Plot settings.
    plotsz = (800, 400)
    fntsz = 12
    fnt = font(fntsz)
    default(
        size=plotsz,
        palette=:Set1_9,
        dpi=300,
        tickfont=fnt,
        guidefont=fnt,
        legendfontsize=fntsz - 2,
        margin=10mm
    )

    pp = plot(xlabel="\$f\$  (GHz)", ylabel="\$|S_{11}|\$  (dB)", legend=:bottomright)
    plot!(pp, f, mag_ref, label="LumpedPort RLC (reference)", lw=2)
    plot!(pp, f, mag_rat, label="RationalImpedance (test)", ls=:dash, lw=2)
    savefig(pp, joinpath(example_dir, "postpro", "rational_impedance.png"))
    display(pp)

    return
end
