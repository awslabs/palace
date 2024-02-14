# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using CSV
using DataFrames
using Measures
using Plots
using PyPlot: matplotlib

"""
    generate_coaxial_data(; num_processors::Integer=1)

Generate the data for the coaxial cable example

# Arguments

  - num_processors - number of processors to use for the simulation
"""
function generate_coaxial_data(; num_processors::Integer=1)
    # Call the solver, discarding the terminal output
    coaxial_dir = @__DIR__
    for sim ∈ ["matched", "open", "short"]
        call_command = Cmd(`palace -np $num_processors coaxial_$sim.json`, dir=coaxial_dir)
        run(call_command)
    end

    # Parse simulation data
    file = joinpath(coaxial_dir, "postpro", "matched", "port-V.csv")
    data_matched = CSV.File(file, header=1) |> DataFrame |> Matrix
    t = data_matched[:, 1]
    data_matched = data_matched[:, 3]
    n_t = size(t, 1)

    file = joinpath(coaxial_dir, "postpro", "open", "port-V.csv")
    data_open = CSV.File(file, header=1) |> DataFrame |> Matrix
    data_open = data_open[:, 3]

    file = joinpath(coaxial_dir, "postpro", "short", "port-V.csv")
    data_short = CSV.File(file, header=1) |> DataFrame |> Matrix
    data_short = data_short[:, 3]

    # Plot settings
    pyplot()
    rcParams = PyPlot.PyDict(matplotlib["rcParams"])
    plotsz = (800, 400)
    fntsz = 12
    fnt = font(fntsz)
    rcParams["mathtext.fontset"] = "stix"
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
