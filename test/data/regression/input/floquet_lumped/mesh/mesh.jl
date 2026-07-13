# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate mesh for the Floquet + lumped port test (from the floquet_lumped directory):
#   julia -e 'include("mesh/mesh.jl"); generate_floquet_lumped_mesh()'

using Gmsh: gmsh

"""
    generate_floquet_lumped_mesh(; filename, ...)

Generate a structured hex mesh for a periodic unit cell with both a Floquet port and a
lumped port, for testing mixed port S-parameter power normalization.

# Geometry

A doubly periodic cell (period `a` in x, `b` in y) with height 2L in z:

  - Floquet port at z = -L (boundary 5)
  - Lumped port resistive gap at z = 0 (boundary 7), oriented in x-direction
  - PEC termination at z = +L (boundary 8)
  - x-periodic boundaries: 1 (x_min), 2 (x_max)
  - y-periodic boundaries: 3 (y_min), 4 (y_max)

The lumped port spans the full cell cross-section (a x b) as a thin horizontal gap at
z = 0, modeling a resistive sheet (like varactors or PIN diodes in a metasurface unit cell).
The Floquet port launches a normally-incident plane wave from below. Power balance
(|S_floquet_reflected|^2 + |S_lumped_absorbed|^2 = 1 for a lossless dielectric) validates
the cross-type normalization bridge.

# Attributes

  - Volume 1: dielectric (below gap, z in [-L, 0])
  - Volume 2: dielectric (above gap, z in [0, +L])
  - Boundary 1,2: x periodic pair
  - Boundary 3,4: y periodic pair
  - Boundary 5: Floquet port (z = -L)
  - Boundary 6: (unused, reserved)
  - Boundary 7: lumped port (z = 0, the full xy cross-section)
  - Boundary 8: PEC (z = +L)

Units: cm (Palace L0 = 1e-2 for cm -> m conversion).
"""
function generate_floquet_lumped_mesh(;
    filename::AbstractString="mesh/floquet_lumped.msh",
    a::Real=1.0,
    b::Real=1.0,
    L::Real=2.0,
    n_x::Integer=8,
    n_y::Integer=8,
    n_z::Integer=20,
    order::Integer=1,
    verbose::Integer=5,
    gui::Bool=false
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    if "floquet_lumped" in gmsh.model.list()
        gmsh.model.setCurrent("floquet_lumped")
        gmsh.model.remove()
    end
    gmsh.model.add("floquet_lumped")

    geo = gmsh.model.geo

    # Simple two-layer box: lower half [-L, 0] and upper half [0, +L]
    # With a shared internal surface at z=0 for the lumped port.

    # Points: 4 corners at each of 3 z-levels (z=-L, z=0, z=+L)
    p = Dict{Tuple{Int, Int, Int}, Int}()
    xs = [0.0, a]
    ys = [0.0, b]
    zs = [-L, 0.0, L]
    for (iz, z) in enumerate(zs), (iy, y) in enumerate(ys), (ix, x) in enumerate(xs)
        p[(ix, iy, iz)] = geo.addPoint(x, y, z)
    end

    # Lines along x (2 per y-level per z-level = 2*2*3 = 12... but only 1 per combo)
    lx = Dict{Tuple{Int, Int, Int}, Int}()
    for iz = 1:3, iy = 1:2
        lx[(1, iy, iz)] = geo.addLine(p[(1, iy, iz)], p[(2, iy, iz)])
    end
    # Lines along y
    ly = Dict{Tuple{Int, Int, Int}, Int}()
    for iz = 1:3, ix = 1:2
        ly[(ix, 1, iz)] = geo.addLine(p[(ix, 1, iz)], p[(ix, 2, iz)])
    end
    # Lines along z (lower half and upper half)
    lz = Dict{Tuple{Int, Int, Int}, Int}()
    for iz = 1:2, iy = 1:2, ix = 1:2
        lz[(ix, iy, iz)] = geo.addLine(p[(ix, iy, iz)], p[(ix, iy, iz + 1)])
    end

    function make_surface(l1, l2, l3, l4)
        cl = geo.addCurveLoop([l1, l2, l3, l4])
        return geo.addPlaneSurface([cl])
    end

    # Z-faces (horizontal): at z=-L, z=0, z=+L
    sz = Dict{Int, Int}()
    for iz = 1:3
        sz[iz] =
            make_surface(lx[(1, 1, iz)], ly[(2, 1, iz)], -lx[(1, 2, iz)], -ly[(1, 1, iz)])
    end

    # X-faces: at x=0 and x=a, for each z-half
    sx = Dict{Tuple{Int, Int}, Int}()
    for ix = 1:2, iz = 1:2
        sx[(ix, iz)] = make_surface(
            ly[(ix, 1, iz)],
            lz[(ix, 2, iz)],
            -ly[(ix, 1, iz + 1)],
            -lz[(ix, 1, iz)]
        )
    end

    # Y-faces: at y=0 and y=b, for each z-half
    sy = Dict{Tuple{Int, Int}, Int}()
    for iy = 1:2, iz = 1:2
        sy[(iy, iz)] = make_surface(
            lx[(1, iy, iz)],
            lz[(2, iy, iz)],
            -lx[(1, iy, iz + 1)],
            -lz[(1, iy, iz)]
        )
    end

    # Volumes: lower (iz=1) and upper (iz=2)
    vols = Dict{Int, Int}()
    for iz = 1:2
        sl = geo.addSurfaceLoop([
            sz[iz],
            sz[iz + 1],
            sx[(1, iz)],
            sx[(2, iz)],
            sy[(1, iz)],
            sy[(2, iz)]
        ])
        vols[iz] = geo.addVolume([sl])
    end

    geo.synchronize()

    # Physical groups
    gmsh.model.addPhysicalGroup(3, [vols[1]], 1, "lower")
    gmsh.model.addPhysicalGroup(3, [vols[2]], 2, "upper")

    gmsh.model.addPhysicalGroup(2, [sx[(1, 1)], sx[(1, 2)]], 1, "x_min")
    gmsh.model.addPhysicalGroup(2, [sx[(2, 1)], sx[(2, 2)]], 2, "x_max")
    gmsh.model.addPhysicalGroup(2, [sy[(1, 1)], sy[(1, 2)]], 3, "y_min")
    gmsh.model.addPhysicalGroup(2, [sy[(2, 1)], sy[(2, 2)]], 4, "y_max")
    gmsh.model.addPhysicalGroup(2, [sz[1]], 5, "floquet_port")
    gmsh.model.addPhysicalGroup(2, [sz[2]], 7, "lumped_port")
    gmsh.model.addPhysicalGroup(2, [sz[3]], 8, "pec")

    # Transfinite constraints for structured mesh
    for iz = 1:3, iy = 1:2
        gmsh.model.mesh.setTransfiniteCurve(lx[(1, iy, iz)], n_x + 1)
    end
    for iz = 1:3, ix = 1:2
        gmsh.model.mesh.setTransfiniteCurve(ly[(ix, 1, iz)], n_y + 1)
    end
    nz_half = max(1, div(n_z, 2))
    for iz = 1:2, iy = 1:2, ix = 1:2
        gmsh.model.mesh.setTransfiniteCurve(lz[(ix, iy, iz)], nz_half + 1)
    end

    for iz = 1:3
        gmsh.model.mesh.setTransfiniteSurface(sz[iz])
        gmsh.model.mesh.setRecombine(2, sz[iz])
    end
    for ix = 1:2, iz = 1:2
        gmsh.model.mesh.setTransfiniteSurface(sx[(ix, iz)])
        gmsh.model.mesh.setRecombine(2, sx[(ix, iz)])
    end
    for iy = 1:2, iz = 1:2
        gmsh.model.mesh.setTransfiniteSurface(sy[(iy, iz)])
        gmsh.model.mesh.setRecombine(2, sy[(iy, iz)])
    end
    for iz = 1:2
        gmsh.model.mesh.setTransfiniteVolume(vols[iz])
    end

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    # Periodic mesh matching
    gmsh.model.mesh.setPeriodic(
        2,
        [sx[(2, 1)], sx[(2, 2)]],
        [sx[(1, 1)], sx[(1, 2)]],
        [1, 0, 0, a, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    )
    gmsh.model.mesh.setPeriodic(
        2,
        [sy[(2, 1)], sy[(2, 2)]],
        [sy[(1, 1)], sy[(1, 2)]],
        [1, 0, 0, 0, 0, 1, 0, b, 0, 0, 1, 0, 0, 0, 0, 1]
    )

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(filename)

    if gui
        gmsh.fltk.run()
    end

    gmsh.finalize()
    return nothing
end
