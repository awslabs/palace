# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate the mesh for the dielectric slab test (from the dielectric_slab directory):
#   julia -e 'include("mesh/mesh.jl"); generate_slab_mesh()'

using Gmsh: gmsh

"""
    generate_slab_mesh(; filename, ...)

Generate a structured hex mesh of a periodic unit cell across which a dielectric slab
extends, for testing the Floquet port S-parameters of the slab at normal incidence against
the Fresnel formula. Without a Floquet wave vector the fields are uniform across the cell,
so a few elements suffice laterally.

# Geometry

A unit cell of period `a` in x and `b` in y, from z = -L to z = L, with a slab of thickness
`h` centered at z = 0.

# Attributes

  - Volume 1: vacuum, Volume 2: slab
  - Boundary 1,2: x periodic pair
  - Boundary 3,4: y periodic pair
  - Boundary 5: port 1 (z = -L), Boundary 6: port 2 (z = +L)

Units: cm (Palace L0 = 1e-2 for cm -> m conversion).
"""
function generate_slab_mesh(;
    filename::AbstractString="mesh/slab.msh",
    a::Real=4.0,
    b::Real=1.0,
    L::Real=4.0,
    h::Real=0.5,
    n_x::Integer=3,
    n_y::Integer=3,
    n_z::Integer=7,
    n_slab_z::Integer=4,
    verbose::Integer=5,
    gui::Bool=false
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    if "slab" in gmsh.model.list()
        gmsh.model.setCurrent("slab")
        gmsh.model.remove()
    end
    gmsh.model.add("slab")

    # Three layers along z (vacuum, slab, vacuum) with shared interfaces.
    zs = [-L, -h / 2, h / 2, L]
    layers = [gmsh.model.occ.addBox(0.0, 0.0, zs[i], a, b, zs[i + 1] - zs[i]) for i = 1:3]
    gmsh.model.occ.fragment([(3, layers[1])], [(3, layers[2]), (3, layers[3])])
    gmsh.model.occ.synchronize()

    # Transfinite constraints for a structured hex mesh: each straight edge is divided
    # according to the axis it runs along, and in z according to its layer.
    for (_, tag) in gmsh.model.getEntities(1)
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(1, tag)
        n = if abs(x1 - x0) > 1.0e-6
            n_x
        elseif abs(y1 - y0) > 1.0e-6
            n_y
        else
            abs(z0 + z1) < h ? n_slab_z : n_z
        end
        gmsh.model.mesh.setTransfiniteCurve(tag, n + 1)
    end
    for (_, tag) in gmsh.model.getEntities(2)
        gmsh.model.mesh.setTransfiniteSurface(tag)
        gmsh.model.mesh.setRecombine(2, tag)
    end
    for (_, tag) in gmsh.model.getEntities(3)
        gmsh.model.mesh.setTransfiniteVolume(tag)
    end

    # Physical groups, from the centers of the volumes and of the outer faces. The faces
    # of each periodic side are sorted by z, so that they pair up in setPeriodic.
    vacuum = Int[]
    slab = Int[]
    for (_, tag) in gmsh.model.getEntities(3)
        z = gmsh.model.occ.getCenterOfMass(3, tag)[3]
        push!(abs(z) < h / 2 ? slab : vacuum, tag)
    end
    sides = [Tuple{Float64, Int}[] for _ = 1:6]
    for (_, tag) in gmsh.model.getEntities(2)
        x, y, z = gmsh.model.occ.getCenterOfMass(2, tag)
        side = if abs(x) < 1.0e-6
            1
        elseif abs(x - a) < 1.0e-6
            2
        elseif abs(y) < 1.0e-6
            3
        elseif abs(y - b) < 1.0e-6
            4
        elseif abs(z + L) < 1.0e-6
            5
        elseif abs(z - L) < 1.0e-6
            6
        else
            0  # Interface between layers
        end
        if side > 0
            push!(sides[side], (z, tag))
        end
    end
    faces = [last.(sort(s)) for s in sides]
    gmsh.model.addPhysicalGroup(3, vacuum, 1, "vacuum")
    gmsh.model.addPhysicalGroup(3, slab, 2, "slab")
    for (i, name) in enumerate(["x_min", "x_max", "y_min", "y_max", "port1", "port2"])
        gmsh.model.addPhysicalGroup(2, faces[i], i, name)
    end

    gmsh.model.mesh.generate(3)

    # Periodic mesh matching: x_max and y_max are translated copies of x_min and y_min.
    gmsh.model.mesh.setPeriodic(
        2,
        faces[2],
        faces[1],
        [1, 0, 0, a, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    )
    gmsh.model.mesh.setPeriodic(
        2,
        faces[4],
        faces[3],
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
