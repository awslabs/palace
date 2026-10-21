# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate mesh for the absorbing boundary waveguide test (from the waveguide_absorbing
# directory):
#   julia -e 'include("mesh/mesh.jl"); generate_waveguide_absorbing_mesh()'

using Gmsh: gmsh

"""
    generate_waveguide_absorbing_mesh(; filename, ...)

Generate a structured hex mesh for a rectangular waveguide terminated by an absorbing
boundary, for testing the reflection of the first- and second-order absorbing boundary
conditions against their closed-form values.

# Geometry

A hollow PEC waveguide of cross-section `a` x `b` and length `L` along z:

  - Wave port at z = 0 (boundary 3), excited in the TE10 mode
  - Absorbing boundary at z = L (boundary 4)
  - PEC side walls (boundary 2)

# Attributes

  - Volume 1: vacuum
  - Boundary 2: PEC side walls
  - Boundary 3: wave port (z = 0)
  - Boundary 4: absorbing boundary (z = L)

Units: mm (Palace L0 = 1e-3 for mm -> m conversion).
"""
function generate_waveguide_absorbing_mesh(;
    filename::AbstractString="mesh/waveguide_absorbing.msh",
    a::Real=20.0,
    b::Real=10.0,
    L::Real=40.0,
    n_x::Integer=4,
    n_y::Integer=2,
    n_z::Integer=8,
    verbose::Integer=5,
    gui::Bool=false
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    if "waveguide_absorbing" in gmsh.model.list()
        gmsh.model.setCurrent("waveguide_absorbing")
        gmsh.model.remove()
    end
    gmsh.model.add("waveguide_absorbing")

    box = gmsh.model.occ.addBox(0.0, 0.0, 0.0, a, b, L)
    gmsh.model.occ.synchronize()

    # Transfinite constraints for a structured hex mesh: each straight edge is divided
    # according to the axis it runs along.
    for (_, tag) in gmsh.model.getEntities(1)
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(1, tag)
        n = abs(x1 - x0) > 1.0e-3 ? n_x : (abs(y1 - y0) > 1.0e-3 ? n_y : n_z)
        gmsh.model.mesh.setTransfiniteCurve(tag, n + 1)
    end
    faces = [tag for (_, tag) in gmsh.model.getBoundary([(3, box)], false, false)]
    for tag in faces
        gmsh.model.mesh.setTransfiniteSurface(tag)
        gmsh.model.mesh.setRecombine(2, tag)
    end
    gmsh.model.mesh.setTransfiniteVolume(box)

    # Physical groups
    port = Int[]
    absorbing = Int[]
    pec = Int[]
    for tag in faces
        z = gmsh.model.occ.getCenterOfMass(2, tag)[3]
        if abs(z) < 1.0e-9
            push!(port, tag)
        elseif abs(z - L) < 1.0e-9
            push!(absorbing, tag)
        else
            push!(pec, tag)
        end
    end
    gmsh.model.addPhysicalGroup(3, [box], 1, "vacuum")
    gmsh.model.addPhysicalGroup(2, pec, 2, "pec")
    gmsh.model.addPhysicalGroup(2, port, 3, "port")
    gmsh.model.addPhysicalGroup(2, absorbing, 4, "absorbing")

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
