# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate compact thin-metal and fabrication-resolved 3D meshes for the
# held-out 20/12 um driven-CPW surface-response validation.

include(joinpath(@__DIR__, "mesh_cpw3d.jl"))

common = (
    w_trace = 20.0,
    w_gap = 12.0,
    w_ground = 40.0,
    h_substrate = 40.0,
    h_vacuum = 40.0,
    lc_far = 12.0,
    lc_transition = 40.0,
    mesh_order = 1,
    extrude_length = 50.0,
    extrude_nz = 2,
)

generate_cpw2d_mesh(;
    common...,
    t_metal = 0.0,
    overetch_depth = 0.0,
    lc_fine = 0.25,
    filename = "cpw3d_surface_maxwell_thin.msh",
)

generate_cpw2d_mesh(;
    common...,
    t_metal = 0.1,
    overetch_depth = 0.05,
    sidewall_angle = 80.0,
    r_top = 0.01,
    r_bot = 0.01,
    lc_fine = 0.01,
    filename = "cpw3d_surface_maxwell_fabricated.msh",
)
