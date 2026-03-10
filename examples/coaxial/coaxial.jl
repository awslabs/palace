# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script processes and visualizes results from the coaxial cable
transient simulation example.

The script runs three Palace simulations with different termination conditions:
1. Matched termination (50Ω load)
2. Open circuit termination
3. Short circuit termination

It then plots the voltage response at the port for all three cases, showing the
effect of reflections from different boundary conditions.

## Prerequisites

This script requires Julia packages. Install them with:

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

## How to run

From the repository root, run:
```bash
julia --project=examples -e 'include("examples/coaxial/coaxial.jl"); generate_coaxial_data(num_processors=4)'
```

This requires `palace` to be a runnable command. If it is not, you can pass the
path to the executable, e.g.
```bash
julia --project=examples -e 'include("examples/coaxial/coaxial.jl"); generate_coaxial_data(palace_exec="build/bin/palace", num_processors=4)'
```

The script will:
1. Run Palace simulations for all three termination cases
2. Parse the voltage data from CSV output files
3. Generate a plot comparing the three cases
4. Save the plot to `postpro/coaxial.png`

## Output

The generated plot shows normalized voltage (V/V₀) vs time for the three
termination conditions, illustrating transmission line reflection behavior.
=#

using CSV
using DataFrames
using Measures
using Plots

"""
    generate_coaxial_data(;
                           palace_exec::String="palace",
                           num_processors::Integer=1
                          )

Generate the data for the coaxial cable example and visualize the results.

# Arguments

  - palace_exec - executable for Palace
  - num_processors - number of processors to use for the simulation
"""
function generate_coaxial_data(; palace_exec="palace", num_processors::Integer=1)
    palace_exec_is_path = occursin(Base.Filesystem.path_separator, palace_exec)
    if palace_exec_is_path
        # Convert palace_exec to absolute path if it's relative
        palace_exec = isabspath(palace_exec) ? palace_exec : abspath(palace_exec)
    end

    # Call the solver, discarding the terminal output
    coaxial_dir = @__DIR__
    for sim ∈ ["matched", "open", "short"]
        base_cmd = `$palace_exec -np $num_processors coaxial_$sim.json`
        call_command = Cmd(base_cmd; dir=coaxial_dir)
        run(call_command)
    end

    # Parse simulation data
    file = joinpath(coaxial_dir, "postpro", "matched", "port-V.csv")
    data_matched = CSV.File(file, header=1) |> DataFrame |> Matrix
    t = data_matched[:, 1]
    data_matched = data_matched[:, 3]

    file = joinpath(coaxial_dir, "postpro", "open", "port-V.csv")
    data_open = CSV.File(file, header=1) |> DataFrame |> Matrix
    data_open = data_open[:, 3]

    file = joinpath(coaxial_dir, "postpro", "short", "port-V.csv")
    data_short = CSV.File(file, header=1) |> DataFrame |> Matrix
    data_short = data_short[:, 3]

    # Plot settings
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

    # Make plots
    xlim = (minimum(t) - 0.1, maximum(t) + 0.1)
    xlbl = "\$t\$  (ns)"
    ylbl = string("\$V\\ /\\ V_0\$")

    pp = plot(xlims=xlim, xlabel=xlbl, ylabel=ylbl, legend=:bottomright)

    lbl = "Open"
    plot!(pp, t, data_open, label=lbl)

    lbl = "Short"
    plot!(pp, t, data_short, label=lbl)

    lbl = "Matched"
    plot!(pp, t, data_matched, label=lbl)

    savefig(pp, joinpath(coaxial_dir, "postpro", "coaxial.png"))
    display(pp)

    return
end
