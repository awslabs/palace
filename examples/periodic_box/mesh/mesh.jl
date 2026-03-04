# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate mesh with:
# julia -e 'include("mesh/mesh.jl"); generate_grating_mesh(filename="mesh/periodic_box.msh")'

using Gmsh: gmsh

"""
    generate_grating_mesh(;
        filename::AbstractString = "periodic_box.msh",
        a::Real = 4.0,
        b::Real = 1.0,
        L::Real = 4.0,
        bar_w::Real = 2.0,
        bar_h::Real = 0.5,
        n_x::Integer = 8,
        n_y::Integer = 2,
        n_z::Integer = 8,
        n_bar_z::Integer = 4,
        order::Integer = 1,
        verbose::Integer = 5,
        gui::Bool = false
    )

Generate a structured hex mesh for a dielectric grating test case.

The geometry is a unit cell with period a in x, b in y, total height 2L in z.
A dielectric bar (centered at z=0, width bar_w, height bar_h) creates a binary
grating that diffracts into multiple orders.

At f=10 GHz (lambda=3cm), a=4cm gives 3 propagating x-orders: (0,0), (±1,0).
The grating bar produces non-trivial diffraction efficiencies.

Layout (cross-section in x-z):
```
  Port 2 (z = L)
  ┌──────────┐ vacuum (eps=1)
  │          │
  │ ████████ │ dielectric bar (eps=4), width=bar_w, centered
  │          │
  └──────────┘ vacuum (eps=1)
  Port 1 (z = -L)
```

Attributes:
  - Volume 1: vacuum
  - Volume 2: dielectric bar
  - Boundary 1,2: x periodic pair
  - Boundary 3,4: y periodic pair
  - Boundary 5: port 1 (z = -L, excitation)
  - Boundary 6: port 2 (z = +L, transmission)

Units: cm (Palace L0 = 1e-2 for cm → m conversion).
"""
function generate_grating_mesh(;
    filename::AbstractString = "periodic_box.msh",
    a::Real = 4.0,
    b::Real = 1.0,
    L::Real = 4.0,
    bar_w::Real = 2.0,
    bar_h::Real = 0.5,
    n_x::Integer = 8,
    n_y::Integer = 2,
    n_z::Integer = 8,
    n_bar_z::Integer = 4,
    order::Integer = 1,
    verbose::Integer = 5,
    gui::Bool = false
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    if "grating" in gmsh.model.list()
        gmsh.model.setCurrent("grating")
        gmsh.model.remove()
    end
    gmsh.model.add("grating")

    geo = gmsh.model.geo

    # The bar is centered in x: from (a-bar_w)/2 to (a+bar_w)/2, and in z: from -bar_h/2
    # to +bar_h/2. We split the domain into 5 z-layers and 3 x-strips.

    x0 = 0.0
    x1 = (a - bar_w) / 2.0  # left edge of bar
    x2 = (a + bar_w) / 2.0  # right edge of bar
    x3 = a

    z0 = -L           # port 1
    z1 = -bar_h / 2.0  # bottom of bar
    z2 = +bar_h / 2.0  # top of bar
    z3 = +L            # port 2

    y0 = 0.0
    y1 = b

    # Number of elements per segment
    nx_left = max(1, round(Int, n_x * (x1 - x0) / a))
    nx_bar = max(1, round(Int, n_x * bar_w / a))
    nx_right = max(1, round(Int, n_x * (x3 - x2) / a))
    nz_bot = max(1, round(Int, n_z * (z1 - z0) / (2 * L)))
    nz_bar = max(1, n_bar_z)
    nz_top = max(1, round(Int, n_z * (z3 - z2) / (2 * L)))

    # Create 9 blocks (3 x-strips × 3 z-layers) using transfinite structured hex.
    # Points: 4 x-positions × 4 z-positions × 2 y-positions = 32 points.
    xs = [x0, x1, x2, x3]
    zs = [z0, z1, z2, z3]
    ys = [y0, y1]
    nxs = [nx_left, nx_bar, nx_right]
    nzs = [nz_bot, nz_bar, nz_top]

    # Create all points
    pts = Dict{Tuple{Int,Int,Int}, Int}()
    for (iz, z) in enumerate(zs)
        for (iy, y) in enumerate(ys)
            for (ix, x) in enumerate(xs)
                pts[(ix, iy, iz)] = geo.addPoint(x, y, z)
            end
        end
    end

    # Create lines in x-direction
    lines_x = Dict{Tuple{Int,Int,Int}, Int}()
    for iz in 1:4, iy in 1:2, ix in 1:3
        lines_x[(ix, iy, iz)] = geo.addLine(pts[(ix, iy, iz)], pts[(ix+1, iy, iz)])
    end

    # Create lines in y-direction
    lines_y = Dict{Tuple{Int,Int,Int}, Int}()
    for iz in 1:4, iy in 1:1, ix in 1:4
        lines_y[(ix, iy, iz)] = geo.addLine(pts[(ix, iy, iz)], pts[(ix, iy+1, iz)])
    end

    # Create lines in z-direction
    lines_z = Dict{Tuple{Int,Int,Int}, Int}()
    for iz in 1:3, iy in 1:2, ix in 1:4
        lines_z[(ix, iy, iz)] = geo.addLine(pts[(ix, iy, iz)], pts[(ix, iy, iz+1)])
    end

    # Helper to create a planar surface from 4 lines (given as oriented line tags)
    function make_surface(l1, l2, l3, l4)
        cl = geo.addCurveLoop([l1, l2, l3, l4])
        return geo.addPlaneSurface([cl])
    end

    # Create all surfaces and volumes for each of the 9 blocks
    volumes = Dict{Tuple{Int,Int}, Int}()
    all_surfaces = []

    for ix in 1:3, iz in 1:3
        # Bottom face (y=y0): lines in x-z plane at y=1
        s_bot = make_surface(
            lines_x[(ix, 1, iz)], lines_z[(ix+1, 1, iz)],
            -lines_x[(ix, 1, iz+1)], -lines_z[(ix, 1, iz)]
        )
        # Top face (y=y1): lines in x-z plane at y=2
        s_top = make_surface(
            lines_x[(ix, 2, iz)], lines_z[(ix+1, 2, iz)],
            -lines_x[(ix, 2, iz+1)], -lines_z[(ix, 2, iz)]
        )
        # Front face (z=z[iz]): lines in x-y plane
        s_front = make_surface(
            lines_x[(ix, 1, iz)], lines_y[(ix+1, 1, iz)],
            -lines_x[(ix, 2, iz)], -lines_y[(ix, 1, iz)]
        )
        # Back face (z=z[iz+1])
        s_back = make_surface(
            lines_x[(ix, 1, iz+1)], lines_y[(ix+1, 1, iz+1)],
            -lines_x[(ix, 2, iz+1)], -lines_y[(ix, 1, iz+1)]
        )
        # Left face (x=x[ix])
        s_left = make_surface(
            lines_y[(ix, 1, iz)], lines_z[(ix, 2, iz)],
            -lines_y[(ix, 1, iz+1)], -lines_z[(ix, 1, iz)]
        )
        # Right face (x=x[ix+1])
        s_right = make_surface(
            lines_y[(ix+1, 1, iz)], lines_z[(ix+1, 2, iz)],
            -lines_y[(ix+1, 1, iz+1)], -lines_z[(ix+1, 1, iz)]
        )

        sl = geo.addSurfaceLoop([s_bot, s_top, s_front, s_back, s_left, s_right])
        volumes[(ix, iz)] = geo.addVolume([sl])
        push!(all_surfaces, (ix, iz, s_bot, s_top, s_front, s_back, s_left, s_right))
    end

    geo.synchronize()

    # Physical groups
    # Vacuum: all blocks except the bar (ix=2, iz=2)
    vac_vols = [volumes[(ix, iz)] for ix in 1:3 for iz in 1:3 if !(ix == 2 && iz == 2)]
    bar_vol = [volumes[(2, 2)]]
    gmsh.model.addPhysicalGroup(3, vac_vols, 1, "vacuum")
    gmsh.model.addPhysicalGroup(3, bar_vol, 2, "dielectric_bar")

    # Boundary groups
    # x_min (x=0) surfaces: ix=1, left faces
    x_min_surfs = [t[7] for t in all_surfaces if t[1] == 1]
    x_max_surfs = [t[8] for t in all_surfaces if t[1] == 3]
    y_min_surfs = [t[3] for t in all_surfaces]  # all front (y=0) faces
    y_max_surfs = [t[4] for t in all_surfaces]  # all back (y=b) faces
    z_min_surfs = [t[5] for t in all_surfaces if t[2] == 1]  # front (z=z0) faces
    z_max_surfs = [t[6] for t in all_surfaces if t[2] == 3]  # back (z=z3) faces

    gmsh.model.addPhysicalGroup(2, x_min_surfs, 1, "x_min")
    gmsh.model.addPhysicalGroup(2, x_max_surfs, 2, "x_max")
    gmsh.model.addPhysicalGroup(2, y_min_surfs, 3, "y_min")
    gmsh.model.addPhysicalGroup(2, y_max_surfs, 4, "y_max")
    gmsh.model.addPhysicalGroup(2, z_min_surfs, 5, "port1")
    gmsh.model.addPhysicalGroup(2, z_max_surfs, 6, "port2")

    # Transfinite constraints
    for iz in 1:4, iy in 1:2
        for (ix, nx) in enumerate(nxs)
            gmsh.model.mesh.setTransfiniteCurve(lines_x[(ix, iy, iz)], nx + 1)
        end
    end
    for iz in 1:4, ix in 1:4
        gmsh.model.mesh.setTransfiniteCurve(lines_y[(1, 1, iz) |> x -> (ix, 1, iz)],
                                             n_y + 1)
    end
    # Simpler: set all y-lines
    for iz in 1:4, ix in 1:4
        gmsh.model.mesh.setTransfiniteCurve(lines_y[(ix, 1, iz)], n_y + 1)
    end
    for iy in 1:2, ix in 1:4
        for (iz, nz) in enumerate(nzs)
            gmsh.model.mesh.setTransfiniteCurve(lines_z[(ix, iy, iz)], nz + 1)
        end
    end

    # Set all surfaces and volumes transfinite
    for (_, _, s1, s2, s3, s4, s5, s6) in all_surfaces
        for s in [s1, s2, s3, s4, s5, s6]
            gmsh.model.mesh.setTransfiniteSurface(s)
            gmsh.model.mesh.setRecombine(2, s)
        end
    end
    for (_, vol) in volumes
        gmsh.model.mesh.setTransfiniteVolume(vol)
    end

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(filename)

    if gui
        gmsh.fltk.run()
    end

    gmsh.finalize()
    return nothing
end
