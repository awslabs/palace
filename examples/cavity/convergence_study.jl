# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
    This script performs a convergence study of the cavity resonator problem

    Loops over mesh refinement level and order
        a) Generates a mesh that corresponds, by calling `mesh.jl`,
        b) Writes a JSON for driving the simulation,
        c) Calls the solver and records the number of DOF,
        d) Extracts from the written CSV files the eigenfrequency.

    Once these are generated, plots DOF^(-1/3) against error.
=#

using DelimitedFiles
using Measures
using Plots
using PyPlot: matplotlib

include(joinpath(@__DIR__, "cavity.jl"))

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

mkrsz = 8
markers = [
    (:circle, mkrsz, stroke(0)),
    (:utriangle, mkrsz, stroke(0)),
    (:square, mkrsz, stroke(0)),
    (:dtriangle, mkrsz, stroke(0)),
    (:star, mkrsz, stroke(0))
]

# Compute the convergence data
num_processors = 32
p_min = 1
p_max = 5
ref_min = 0
ref_max = 3
for mesh_type ∈ [0, 1, 2]
    if mesh_type == 0
        mesh_name = "Tetrahedra"
    elseif mesh_type == 1
        mesh_name = "Prism"
    elseif mesh_type == 2
        mesh_name = "Hexahedra"
    end
    println("$mesh_name:")

    # Run simulations
    dof, f_TM_010_relative_error, f_TE_111_relative_error =
        generate_cavity_convergence_data(
            p_min=p_min,
            p_max=p_max,
            ref_min=ref_min,
            ref_max=ref_max,
            mesh_type=mesh_type,
            num_processors=num_processors
        )

    # Plot the convergence
    xlbl = "\$DOF^{-1/3}\$"
    ylbl = "Relative error, $mesh_name"
    pp = plot(xlabel=xlbl, ylabel=ylbl, legend=:bottomright)
    for p ∈ p_min:p_max
        plot!(
            pp,
            dof[p] .^ (-1 / 3),
            f_TM_010_relative_error[p],
            label=string("\$f^{TM}_{010}, p = ", p, "\$"),
            linestyle=:solid,
            markers=markers[p],
            color=p
        )
        plot!(
            pp,
            dof[p] .^ (-1 / 3),
            f_TE_111_relative_error[p],
            label=string("\$f^{TE}_{111}, p = ", p, "\$"),
            linestyle=:dash,
            markers=markers[p],
            color=p
        )
    end
    plot!(pp, xaxis=:log, yaxis=:log)

    # Compute the rate from the final entries in the relative error
    # Let that e ~ C * h^k, where h ~ DOF^(-1/3), then log and compute the slopes between
    # points
    Δlogh = map(x -> log.(x[2:end] .^ (-1 / 3)) - log.(x[1:(end - 1)] .^ (-1 / 3)), dof)
    Δlogf_TM_010 = map(x -> log.(x[2:end]) - log.(x[1:(end - 1)]), f_TM_010_relative_error)
    Δlogf_TE_111 = map(x -> log.(x[2:end]) - log.(x[1:(end - 1)]), f_TE_111_relative_error)

    k_f_TM_010 = map((x, y) -> x ./ y, Δlogf_TM_010, Δlogh)
    k_f_TE_111 = map((x, y) -> x ./ y, Δlogf_TE_111, Δlogh)

    println("k_f_TM_010 =", map(x -> round.(x, digits=2), k_f_TM_010))
    println("k_f_TE_111 =", map(x -> round.(x, digits=2), k_f_TE_111))

    cavity_dir = @__DIR__
    output_dir = joinpath(cavity_dir, "postpro", "convergence")
    lmesh_name = lowercase(mesh_name)

    savefig(pp, joinpath(output_dir, string("cavity_error_", lmesh_name, ".png")))
    display(pp)

    writedlm(joinpath(output_dir, string("dof_", lmesh_name, ".csv")), dof)
    writedlm(
        joinpath(output_dir, string("f_TM_010_error_", lmesh_name, ".csv")),
        f_TM_010_relative_error
    )
    writedlm(
        joinpath(output_dir, string("f_TE_111_error_", lmesh_name, ".csv")),
        f_TE_111_relative_error
    )
end
