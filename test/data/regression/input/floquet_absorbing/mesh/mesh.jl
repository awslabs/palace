# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate meshes for the Floquet absorbing boundary tests (from the floquet_absorbing
# directory):
#   julia -e 'include("mesh/mesh.jl"); generate_floquet_absorbing_mesh();
#             generate_floquet_absorbing_mesh(filename="mesh/floquet_absorbing_scatterer.msh",
#                                             scatterer=true)'

using Gmsh: gmsh

"""
    generate_floquet_absorbing_mesh(; filename, ...)

Generate a structured hex mesh of a box with one boundary attribute per face, for testing
the second-order absorbing boundary condition with Floquet periodic boundary conditions.

# Geometry

A vacuum box of `a` x `b` x `L` along x, y, and z. The faces normal to x and y are either
periodic pairs or PEC, the face z = 0 is a Floquet port or PEC, and the face z = L is the
absorbing boundary. With `scatterer = true`, the block `0 < x < a / 2`, `L / 4 < z < L / 2`
is a separate volume, for a dielectric scatterer, so that the periodic part of the field
varies in the cell. The higher Floquet orders it scatters decay before they reach the
absorbing boundary, where the periodic field is nearly uniform.

# Attributes

  - Volume 1: vacuum
  - Volume 2: scatterer (with `scatterer = true`)
  - Boundary 1, 2: x = 0, x = a
  - Boundary 3, 4: y = 0, y = b
  - Boundary 5: z = 0
  - Boundary 6: z = L (absorbing boundary)

Units: cm (Palace L0 = 1e-2 for cm -> m conversion).
"""
function generate_floquet_absorbing_mesh(;
    filename::AbstractString="mesh/floquet_absorbing.msh",
    a::Real=2.0,
    b::Real=2.0,
    L::Real=4.0,
    n_x::Integer=4,
    n_y::Integer=6,
    n_z::Integer=8,
    scatterer::Bool=false,
    verbose::Integer=5,
    gui::Bool=false
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    if "floquet_absorbing" in gmsh.model.list()
        gmsh.model.setCurrent("floquet_absorbing")
        gmsh.model.remove()
    end
    gmsh.model.add("floquet_absorbing")

    # Boxes of a structured block decomposition, conformal after fragmentation. Every box
    # spans the full period in y.
    xs = scatterer ? [0.0, a / 2, a] : [0.0, a]
    zs = scatterer ? [0.0, L / 4, L / 2, L] : [0.0, L]
    boxes = Int[]
    for i = 1:(length(xs) - 1), k = 1:(length(zs) - 1)
        dx, dz = xs[i + 1] - xs[i], zs[k + 1] - zs[k]
        push!(boxes, gmsh.model.occ.addBox(xs[i], 0.0, zs[k], dx, b, dz))
    end
    if length(boxes) > 1
        gmsh.model.occ.fragment([(3, tag) for tag in boxes], [])
    end
    gmsh.model.occ.synchronize()

    # Transfinite constraints for a structured hex mesh: each straight edge is divided
    # according to the axis it runs along and its length. The periodic faces have matching
    # meshes.
    for (_, tag) in gmsh.model.getEntities(1)
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(1, tag)
        n = if abs(x1 - x0) > 1.0e-3
            round(Int, n_x * abs(x1 - x0) / a)
        elseif abs(y1 - y0) > 1.0e-3
            n_y
        else
            round(Int, n_z * abs(z1 - z0) / L)
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

    # Physical groups: volumes by position and one boundary attribute per face of the box.
    vacuum = Int[]
    block = Int[]
    for (_, tag) in gmsh.model.getEntities(3)
        x, y, z = gmsh.model.occ.getCenterOfMass(3, tag)
        push!(scatterer && x < a / 2 && L / 4 < z < L / 2 ? block : vacuum, tag)
    end
    gmsh.model.addPhysicalGroup(3, vacuum, 1, "vacuum")
    if scatterer
        gmsh.model.addPhysicalGroup(3, block, 2, "scatterer")
    end
    faces = Dict(attr => Int[] for attr = 1:6)
    for (_, tag) in gmsh.model.getEntities(2)
        x, y, z = gmsh.model.occ.getCenterOfMass(2, tag)
        if abs(x) < 1.0e-9
            push!(faces[1], tag)
        elseif abs(x - a) < 1.0e-9
            push!(faces[2], tag)
        elseif abs(y) < 1.0e-9
            push!(faces[3], tag)
        elseif abs(y - b) < 1.0e-9
            push!(faces[4], tag)
        elseif abs(z) < 1.0e-9
            push!(faces[5], tag)
        elseif abs(z - L) < 1.0e-9
            push!(faces[6], tag)
        end
    end
    for attr = 1:6
        gmsh.model.addPhysicalGroup(2, faces[attr], attr)
    end

    gmsh.model.mesh.generate(3)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(filename)

    if gui
        gmsh.fltk.run()
    end

    gmsh.finalize()
    return nothing
end
