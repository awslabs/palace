# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate mesh from examples/pml_oblique:
#   julia --project=.. -e 'include("mesh/mesh.jl"); generate_pml_oblique_mesh()'

using Gmsh: gmsh

"""
    generate_pml_oblique_mesh(; filename, a, b, L, d, n_x, n_y, n_z, n_pml_z,
                              order, verbose, gui)

Generate a structured periodic vacuum cell for oblique-incidence PML tests.

The domain is a rectangular unit cell, periodic in x/y. A Floquet port on z = 0
injects a plane wave into a physical vacuum region followed by a +z PML slab and
a PEC termination.

Attributes:

  - Volume 1: physical vacuum
  - Volume 2: +z PML slab
  - Boundary 3,4: x periodic pair
  - Boundary 5,6: y periodic pair
  - Boundary 7: Floquet port, z = 0
  - Boundary 8: outer PML termination, z = L + d

Units are meters (Palace L0 = 1).
"""
function generate_pml_oblique_mesh(;
    filename::AbstractString="mesh/pml_oblique.msh",
    a::Real=0.04,
    b::Real=0.04,
    L::Real=0.20,
    d::Real=0.30,
    n_x::Integer=4,
    n_y::Integer=4,
    n_z::Integer=16,
    n_pml_z::Integer=24,
    order::Integer=1,
    verbose::Integer=4,
    gui::Bool=false
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    if "pml_oblique" in gmsh.model.list()
        gmsh.model.setCurrent("pml_oblique")
        gmsh.model.remove()
    end
    gmsh.model.add("pml_oblique")

    geo = gmsh.model.geo
    xs = [0.0, a]
    ys = [0.0, b]
    zs = [0.0, L, L + d]

    pts = Dict{Tuple{Int, Int, Int}, Int}()
    for (iz, z) in enumerate(zs), (iy, y) in enumerate(ys), (ix, x) in enumerate(xs)
        pts[(ix, iy, iz)] = geo.addPoint(x, y, z)
    end

    lines_x = Dict{Tuple{Int, Int, Int}, Int}()
    for iz = 1:3, iy = 1:2
        lines_x[(1, iy, iz)] = geo.addLine(pts[(1, iy, iz)], pts[(2, iy, iz)])
    end
    lines_y = Dict{Tuple{Int, Int, Int}, Int}()
    for iz = 1:3, ix = 1:2
        lines_y[(ix, 1, iz)] = geo.addLine(pts[(ix, 1, iz)], pts[(ix, 2, iz)])
    end
    lines_z = Dict{Tuple{Int, Int, Int}, Int}()
    for iz = 1:2, iy = 1:2, ix = 1:2
        lines_z[(ix, iy, iz)] = geo.addLine(pts[(ix, iy, iz)], pts[(ix, iy, iz + 1)])
    end

    function make_surface(l1, l2, l3, l4)
        cl = geo.addCurveLoop([l1, l2, l3, l4])
        return geo.addPlaneSurface([cl])
    end

    surfs_x = Dict{Tuple{Int, Int}, Int}()
    for ix = 1:2, iz = 1:2
        surfs_x[(ix, iz)] = make_surface(
            lines_y[(ix, 1, iz)],
            lines_z[(ix, 2, iz)],
            -lines_y[(ix, 1, iz + 1)],
            -lines_z[(ix, 1, iz)]
        )
    end

    surfs_y = Dict{Tuple{Int, Int}, Int}()
    for iy = 1:2, iz = 1:2
        surfs_y[(iy, iz)] = make_surface(
            lines_x[(1, iy, iz)],
            lines_z[(2, iy, iz)],
            -lines_x[(1, iy, iz + 1)],
            -lines_z[(1, iy, iz)]
        )
    end

    surfs_z = Dict{Int, Int}()
    for iz = 1:3
        surfs_z[iz] = make_surface(
            lines_x[(1, 1, iz)],
            lines_y[(2, 1, iz)],
            -lines_x[(1, 2, iz)],
            -lines_y[(1, 1, iz)]
        )
    end

    volumes = Dict{Int, Int}()
    for iz = 1:2
        sl = geo.addSurfaceLoop([
            -surfs_x[(1, iz)],
            surfs_x[(2, iz)],
            surfs_y[(1, iz)],
            -surfs_y[(2, iz)],
            -surfs_z[iz],
            surfs_z[iz + 1]
        ])
        volumes[iz] = geo.addVolume([sl])
    end

    geo.synchronize()

    gmsh.model.addPhysicalGroup(3, [volumes[1]], 1, "vacuum")
    gmsh.model.addPhysicalGroup(3, [volumes[2]], 2, "pml_z_plus")
    gmsh.model.addPhysicalGroup(2, [surfs_x[(1, iz)] for iz = 1:2], 3, "x_min")
    gmsh.model.addPhysicalGroup(2, [surfs_x[(2, iz)] for iz = 1:2], 4, "x_max")
    gmsh.model.addPhysicalGroup(2, [surfs_y[(1, iz)] for iz = 1:2], 5, "y_min")
    gmsh.model.addPhysicalGroup(2, [surfs_y[(2, iz)] for iz = 1:2], 6, "y_max")
    gmsh.model.addPhysicalGroup(2, [surfs_z[1]], 7, "floquet_port")
    gmsh.model.addPhysicalGroup(2, [surfs_z[3]], 8, "outer_pec")

    for iz = 1:3, iy = 1:2
        gmsh.model.mesh.setTransfiniteCurve(lines_x[(1, iy, iz)], n_x + 1)
    end
    for iz = 1:3, ix = 1:2
        gmsh.model.mesh.setTransfiniteCurve(lines_y[(ix, 1, iz)], n_y + 1)
    end
    for iy = 1:2, ix = 1:2
        gmsh.model.mesh.setTransfiniteCurve(lines_z[(ix, iy, 1)], n_z + 1)
        gmsh.model.mesh.setTransfiniteCurve(lines_z[(ix, iy, 2)], n_pml_z + 1)
    end

    for dct in (surfs_x, surfs_y)
        for (_, s) in dct
            gmsh.model.mesh.setTransfiniteSurface(s)
            gmsh.model.mesh.setRecombine(2, s)
        end
    end
    for (_, s) in surfs_z
        gmsh.model.mesh.setTransfiniteSurface(s)
        gmsh.model.mesh.setRecombine(2, s)
    end
    for (_, vol) in volumes
        gmsh.model.mesh.setTransfiniteVolume(vol)
    end

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    gmsh.model.mesh.setPeriodic(
        2,
        [surfs_x[(2, iz)] for iz = 1:2],
        [surfs_x[(1, iz)] for iz = 1:2],
        [1, 0, 0, a, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    )
    gmsh.model.mesh.setPeriodic(
        2,
        [surfs_y[(2, iz)] for iz = 1:2],
        [surfs_y[(1, iz)] for iz = 1:2],
        [1, 0, 0, 0, 0, 1, 0, b, 0, 0, 1, 0, 0, 0, 0, 1]
    )

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(filename)

    println("Generated $(filename)")
    println("  Periods: a = $(a) m, b = $(b) m")
    println("  Physical length: $(L) m, PML thickness: $(d) m")

    if gui
        gmsh.fltk.run()
    end
    gmsh.finalize()
    return nothing
end
