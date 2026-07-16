# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate matched thin-metal and fabrication-resolved CPW meshes for the
# surface-participation calibration study.

include(joinpath(@__DIR__, "mesh_surface_calibration.jl"))

common = (
    w_trace = 10.0,
    w_gap = 6.0,
    w_ground = 300.0,
    h_substrate = 525.0,
    h_vacuum = 500.0,
    lc_far = 40.0,
    mesh_order = 2,
)

generate_cpw2d_mesh(;
    common...,
    t_metal_L1 = 0.0,
    overetch_L1 = 0.0,
    lc_fine = 0.03,
    filename = "cpw2d_surface_calibration_thin.msh",
)

generate_cpw2d_mesh(;
    common...,
    t_metal_L1 = 0.1,
    overetch_L1 = 0.05,
    sidewall_angle = 80.0,
    r_top = 0.01,
    r_bot = 0.01,
    lc_fine = 0.01,
    filename = "cpw2d_surface_calibration_fabricated.msh",
)

validation = merge(common, (w_trace = 20.0, w_gap = 12.0, lc_far = 80.0))

generate_cpw2d_mesh(;
    validation...,
    t_metal_L1 = 0.0,
    overetch_L1 = 0.0,
    lc_fine = 0.03,
    filename = "cpw2d_surface_validation_thin.msh",
)

generate_cpw2d_mesh(;
    validation...,
    t_metal_L1 = 0.1,
    overetch_L1 = 0.05,
    sidewall_angle = 80.0,
    r_top = 0.01,
    r_bot = 0.01,
    lc_fine = 0.01,
    filename = "cpw2d_surface_validation_fabricated.msh",
)
