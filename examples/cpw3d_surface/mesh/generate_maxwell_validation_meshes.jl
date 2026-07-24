# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate compact thin-metal and fabrication-resolved 3D meshes for the
# held-out 20/12 um Maxwell CPW surface-response validation.

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
)

for length in (50.0, 200.0)
    suffix = length == 50.0 ? "" : "_L$(Int(length))"
    extrude_nz = Int(length / 25.0)

    generate_cpw2d_mesh(;
        common...,
        t_metal = 0.0,
        overetch_depth = 0.0,
        lc_fine = 0.25,
        extrude_length = length,
        extrude_nz = extrude_nz,
        filename = "cpw3d_surface_maxwell_thin$(suffix).msh",
    )

    generate_cpw2d_mesh(;
        common...,
        t_metal = 0.1,
        overetch_depth = 0.05,
        sidewall_angle = 80.0,
        r_top = 0.01,
        r_bot = 0.01,
        lc_fine = 0.01,
        extrude_length = length,
        extrude_nz = extrude_nz,
        filename = "cpw3d_surface_maxwell_fabricated$(suffix).msh",
    )
end
