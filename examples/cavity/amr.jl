# SPDX-License-Identifier: Apache-2.0

#=
    This script performs an eigenmode cavity simulation using adaptive mesh refinement, and explores the results.
=#

# Paths
include(joinpath(@__DIR__, "cavity.jl"))
cavity_dir = @__DIR__
output_dirname = output_dir = joinpath(cavity_dir, "postpro", "amr")

# Parameters
num_processors = 4
amr_tol = 1E-2
amr_max_its = 10 # something larger than what is required to hit the error tolerance
amr_max_size = 2E6 # something larger than what is required to hit the error tolerance
amr_update_fraction = 0.7
amr_nonconformal = true

# Run simulations
generate_cavity_amr_data(
    p=2, # use second order (can modify)
    mesh_type=1, # use prism
    num_processors=num_processors,
    amr_max_its=amr_max_its,
    amr_tol=amr_tol,
    amr_max_size=amr_max_size,
    amr_nonconformal=amr_nonconformal,
    amr_update_fraction=amr_update_fraction,
    dirname=output_dir
)
