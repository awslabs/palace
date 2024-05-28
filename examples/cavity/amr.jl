# SPDX-License-Identifier: Apache-2.0

#=
    This script performs an eigenmode cavity simulation using adaptive mesh refinement, and explores the results.
=#

# Paths
include(joinpath(@__DIR__, "cavity.jl"))
cavity_dir = @__DIR__
output_dir_conformal = joinpath(cavity_dir, "postpro", "amr", "conformal")
output_dir_nonconformal = joinpath(cavity_dir, "postpro", "amr", "nonconformal")

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
