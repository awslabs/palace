# SPDX-License-Identifier: Apache-2.0

#=
    This script performs an eigenmode cavity simulation using adaptive mesh refinement, and explores the results.
=#

using DelimitedFiles
using Measures
using Plots
using PyPlot: matplotlib

# Paths
include(joinpath(@__DIR__, "cavity.jl"))
cavity_dir = @__DIR__
output_dir_conformal = joinpath(cavity_dir, "postpro", "amr", "conformal")
output_dir_nonconformal = joinpath(cavity_dir, "postpro", "amr", "nonconformal")
amr_error_plot_file = joinpath(cavity_dir, "postpro", "amr", "cavity_amr_error.png")

# Parameters
num_processors = 8
amr_tol = 1E-2 # something not too aggressive for the sake of speed
amr_max_its = 10 # something larger than what is required to hit the error tolerance
amr_max_size = 2E6 # something larger than what is required to hit the error tolerance
amr_update_fraction = 0.7

p = 2
radius = 2.74
aspect_ratio = 1 / sqrt(2)

# Run simulations
for (amr_nonconformal, output_dir) in
    [(true, output_dir_nonconformal), (false, output_dir_conformal)]
    generate_cavity_amr_data(
        radius=radius,
        aspect_ratio=aspect_ratio,
        p=p, # use second order (can modify)
        mesh_type=0, # use tets (cannot use wedge elements for conformal)
        num_processors=num_processors,
        amr_max_its=amr_max_its,
        amr_tol=amr_tol,
        amr_max_size=amr_max_size,
        amr_nonconformal=amr_nonconformal,
        amr_update_fraction=amr_update_fraction,
        dirname=output_dir
    )
end

# Compute the exact solution for reference
~, f_TM_010_true = frequency_transverse(
    0,
    1,
    0;
    ϵᵣ   = 2.08,
    μᵣ   = 1.0,
    a_cm = radius,
    d_cm = aspect_ratio * 2 * radius
)
f_TE_111_true, ~ = frequency_transverse(
    1,
    1,
    1;
    ϵᵣ   = 2.08,
    μᵣ   = 1.0,
    a_cm = radius,
    d_cm = aspect_ratio * 2 * radius
)

# Parse the results
# go through iterations to look for DOFs and compare to analytical answer
# Loop through each subdirectory in the "amr" directory and process JSON and CSV files

dof_conformal = Vector{Int}()
f_TM_010_conformal = Vector{Float64}()
f_TE_111_conformal = Vector{Float64}()
err_norm_conformal = Vector{Float64}()

dof_nonconformal = Vector{Int}()
f_TM_010_nonconformal = Vector{Float64}()
f_TE_111_nonconformal = Vector{Float64}()
err_norm_nonconformal = Vector{Float64}()

for (dof, f_TM_010, f_TE_111, err_norm, output_dir) in [
    (
        dof_conformal,
        f_TM_010_conformal,
        f_TE_111_conformal,
        err_norm_conformal,
        output_dir_conformal
    ),
    (
        dof_nonconformal,
        f_TM_010_nonconformal,
        f_TE_111_nonconformal,
        err_norm_nonconformal,
        output_dir_nonconformal
    )
]
    subdirectories = [filter(isdir, readdir(output_dir, join=true)); output_dir]

    for subdir in subdirectories
        # Read degrees of freedom from log
        json_files = filter(f -> occursin(r"\.json$", f), readdir(subdir, join=true))
        for json_file in json_files
            json_data = JSON.parsefile(json_file)
            push!(dof, json_data["Problem"]["DegreesOfFreedom"])
        end
        # Read frequencies from csv, and process
        eig_files = filter(f -> occursin("eig.csv", f), readdir(subdir, join=true))
        for eig_file in eig_files
            eig_data = CSV.read(eig_file, DataFrame)
            rename!(eig_data, strip.(names(eig_data)))
            frequencies = Float64[]  # Temporary array to store frequencies for this file
            for row in eachrow(eig_data)
                push!(frequencies, row["Re{f} (GHz)"])
            end
            # The first dominant frequency should be the magnetic, and the second the
            # electric, but just in case we search for the closest
            push!(f_TM_010, frequencies[argmin(abs.(frequencies .- f_TM_010_true))])
            push!(f_TE_111, frequencies[argmin(abs.(frequencies .- f_TE_111_true))])
        end
        # Read error markers from csv
        error_files =
            filter(f -> occursin("error-indicators.csv", f), readdir(subdir, join=true))
        for error_file in error_files
            error_data = CSV.read(error_file, DataFrame)
            rename!(error_data, strip.(names(error_data)))
            for row in eachrow(error_data)
                push!(err_norm, row["Norm"])
            end
        end
    end
end

#=
    Generic plot settings
=#
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

#=
    AME error plot
=#

xlbl = "Iteration number"
ylbl = "AMR Error Norm"
pp = plot(xlabel=xlbl, ylabel=ylbl, legend=:topright)

# Plot error norms vs iteration number
plot!(
    pp,
    1:length(err_norm_conformal),
    err_norm_conformal,
    label=string("p=", p, ", conformal"),
    markers=:circle,
    color=:green
)
plot!(
    pp,
    1:length(err_norm_nonconformal),
    err_norm_nonconformal,
    label=string("p=", p, ", nonconformal"),
    markers=:square,
    color=:purple
)

hline!(pp, [amr_tol], label="AMR Tolerance", linestyle=:dash, color=:black, alpha=0.5)
savefig(pp, amr_error_plot_file)
