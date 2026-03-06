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
        bar_d::Real = 0.5,
        bar_h::Real = 0.5,
        n_x::Integer = 8,
        n_y::Integer = 8,
        n_z::Integer = 8,
        n_bar_z::Integer = 4,
        order::Integer = 1,
        verbose::Integer = 5,
        gui::Bool = false
    )

Generate a structured hex mesh for a 3D dielectric grating test case.

The geometry is a unit cell with period a in x, b in y, total height 2L in z.
A dielectric bar (centered at the origin, width bar_w in x, depth bar_d in y,
height bar_h in z) creates a 3D grating that diffracts into multiple orders.

The domain is split into 27 blocks (3 × 3 × 3) with the bar occupying the center
block. Each block face is created exactly once and shared between adjacent volumes
via signed surface loop tags, ensuring correct boundary element export for MPI.

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
    bar_d::Real = 0.5,
    bar_h::Real = 0.5,
    n_x::Integer = 8,
    n_y::Integer = 8,
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

    # Block boundaries: 4 positions in each direction (3 strips each)
    xs = [0.0, (a - bar_w) / 2.0, (a + bar_w) / 2.0, a]
    ys = [0.0, (b - bar_d) / 2.0, (b + bar_d) / 2.0, b]
    zs = [-L, -bar_h / 2.0, bar_h / 2.0, L]

    # Element counts per segment
    nxs = [max(1, round(Int, n_x * (xs[2] - xs[1]) / a)),
           max(1, round(Int, n_x * (xs[3] - xs[2]) / a)),
           max(1, round(Int, n_x * (xs[4] - xs[3]) / a))]
    nys = [max(1, round(Int, n_y * (ys[2] - ys[1]) / b)),
           max(1, round(Int, n_y * (ys[3] - ys[2]) / b)),
           max(1, round(Int, n_y * (ys[4] - ys[3]) / b))]
    nzs = [max(1, round(Int, n_z * (zs[2] - zs[1]) / (2 * L))),
           max(1, n_bar_z),
           max(1, round(Int, n_z * (zs[4] - zs[3]) / (2 * L)))]

    # --- Points: 4 × 4 × 4 = 64 ---
    pts = Dict{Tuple{Int,Int,Int},Int}()
    for (iz, z) in enumerate(zs), (iy, y) in enumerate(ys), (ix, x) in enumerate(xs)
        pts[(ix, iy, iz)] = geo.addPoint(x, y, z)
    end

    # --- Lines: 48 + 48 + 48 = 144 ---
    lines_x = Dict{Tuple{Int,Int,Int},Int}()
    for iz in 1:4, iy in 1:4, ix in 1:3
        lines_x[(ix, iy, iz)] = geo.addLine(pts[(ix, iy, iz)], pts[(ix + 1, iy, iz)])
    end

    lines_y = Dict{Tuple{Int,Int,Int},Int}()
    for iz in 1:4, iy in 1:3, ix in 1:4
        lines_y[(ix, iy, iz)] = geo.addLine(pts[(ix, iy, iz)], pts[(ix, iy + 1, iz)])
    end

    lines_z = Dict{Tuple{Int,Int,Int},Int}()
    for iz in 1:3, iy in 1:4, ix in 1:4
        lines_z[(ix, iy, iz)] = geo.addLine(pts[(ix, iy, iz)], pts[(ix, iy, iz + 1)])
    end

    function make_surface(l1, l2, l3, l4)
        cl = geo.addCurveLoop([l1, l2, l3, l4])
        return geo.addPlaneSurface([cl])
    end

    # --- 108 unique surfaces ---
    # Surface normals (right-hand rule on curve loop):
    #   surfs_x: +x    surfs_y: -y    surfs_z: +z

    # X-faces: 4 × 3 × 3 = 36
    surfs_x = Dict{Tuple{Int,Int,Int},Int}()
    for ix in 1:4, iy in 1:3, iz in 1:3
        surfs_x[(ix, iy, iz)] = make_surface(
            lines_y[(ix, iy, iz)], lines_z[(ix, iy + 1, iz)],
            -lines_y[(ix, iy, iz + 1)], -lines_z[(ix, iy, iz)]
        )
    end

    # Y-faces: 3 × 4 × 3 = 36
    surfs_y = Dict{Tuple{Int,Int,Int},Int}()
    for ix in 1:3, iy in 1:4, iz in 1:3
        surfs_y[(ix, iy, iz)] = make_surface(
            lines_x[(ix, iy, iz)], lines_z[(ix + 1, iy, iz)],
            -lines_x[(ix, iy, iz + 1)], -lines_z[(ix, iy, iz)]
        )
    end

    # Z-faces: 3 × 3 × 4 = 36
    surfs_z = Dict{Tuple{Int,Int,Int},Int}()
    for ix in 1:3, iy in 1:3, iz in 1:4
        surfs_z[(ix, iy, iz)] = make_surface(
            lines_x[(ix, iy, iz)], lines_y[(ix + 1, iy, iz)],
            -lines_x[(ix, iy + 1, iz)], -lines_y[(ix, iy, iz)]
        )
    end

    # --- 27 volumes from shared surfaces ---
    # Outward normals for volume (ix, iy, iz):
    #   left  (x[ix]):    need -x → negate surfs_x (natural +x)
    #   right (x[ix+1]):  need +x → keep
    #   front (y[iy]):    need -y → surfs_y is -y, keep
    #   back  (y[iy+1]):  need +y → negate surfs_y
    #   bottom (z[iz]):   need -z → negate surfs_z (natural +z)
    #   top    (z[iz+1]): need +z → keep
    volumes = Dict{Tuple{Int,Int,Int},Int}()
    for ix in 1:3, iy in 1:3, iz in 1:3
        sl = geo.addSurfaceLoop([
            -surfs_x[(ix, iy, iz)],
             surfs_x[(ix + 1, iy, iz)],
             surfs_y[(ix, iy, iz)],
            -surfs_y[(ix, iy + 1, iz)],
            -surfs_z[(ix, iy, iz)],
             surfs_z[(ix, iy, iz + 1)]
        ])
        volumes[(ix, iy, iz)] = geo.addVolume([sl])
    end

    geo.synchronize()

    # --- Physical groups ---
    vac_vols = [volumes[(ix, iy, iz)]
                for ix in 1:3 for iy in 1:3 for iz in 1:3
                if !(ix == 2 && iy == 2 && iz == 2)]
    gmsh.model.addPhysicalGroup(3, vac_vols, 1, "vacuum")
    gmsh.model.addPhysicalGroup(3, [volumes[(2, 2, 2)]], 2, "dielectric_bar")

    gmsh.model.addPhysicalGroup(2,
        [surfs_x[(1, iy, iz)] for iy in 1:3 for iz in 1:3], 1, "x_min")
    gmsh.model.addPhysicalGroup(2,
        [surfs_x[(4, iy, iz)] for iy in 1:3 for iz in 1:3], 2, "x_max")
    gmsh.model.addPhysicalGroup(2,
        [surfs_y[(ix, 1, iz)] for ix in 1:3 for iz in 1:3], 3, "y_min")
    gmsh.model.addPhysicalGroup(2,
        [surfs_y[(ix, 4, iz)] for ix in 1:3 for iz in 1:3], 4, "y_max")
    gmsh.model.addPhysicalGroup(2,
        [surfs_z[(ix, iy, 1)] for ix in 1:3 for iy in 1:3], 5, "port1")
    gmsh.model.addPhysicalGroup(2,
        [surfs_z[(ix, iy, 4)] for ix in 1:3 for iy in 1:3], 6, "port2")

    # --- Transfinite constraints ---
    for iz in 1:4, iy in 1:4
        for (ix, nx) in enumerate(nxs)
            gmsh.model.mesh.setTransfiniteCurve(lines_x[(ix, iy, iz)], nx + 1)
        end
    end
    for iz in 1:4, ix in 1:4
        for (iy, ny) in enumerate(nys)
            gmsh.model.mesh.setTransfiniteCurve(lines_y[(ix, iy, iz)], ny + 1)
        end
    end
    for iy in 1:4, ix in 1:4
        for (iz, nz) in enumerate(nzs)
            gmsh.model.mesh.setTransfiniteCurve(lines_z[(ix, iy, iz)], nz + 1)
        end
    end

    for d in (surfs_x, surfs_y, surfs_z)
        for (_, s) in d
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
