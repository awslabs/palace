# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate fabrication-resolved and thin-metal 3D extrusions of the held-out
# 20/12 um CPW used by the 2D surface-calibration validation study.

include(joinpath(@__DIR__, "mesh_cpw3d.jl"))

common = (
    w_trace = 20.0,
    w_gap = 12.0,
    w_ground = 300.0,
    h_substrate = 525.0,
    h_vacuum = 500.0,
    lc_far = 80.0,
    mesh_order = 2,
    extrude_nz = 4,
)

for length in (50.0, 200.0)
    generate_cpw2d_mesh(;
        common...,
        t_metal = 0.0,
        overetch_depth = 0.0,
        lc_fine = 0.03,
        extrude_length = length,
        filename = "cpw3d_surface_validation_thin_L$(Int(length)).msh",
    )
end

generate_cpw2d_mesh(;
    common...,
    t_metal = 0.1,
    overetch_depth = 0.05,
    sidewall_angle = 80.0,
    r_top = 0.01,
    r_bot = 0.01,
    lc_fine = 0.01,
    extrude_length = 50.0,
    filename = "cpw3d_surface_validation_fabricated_L50.msh",
)
