# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate the mesh of the silicon carbide half-space (from the dielectric_grating example
# directory):
#   julia -e 'include("mesh/slab.jl"); generate_slab_mesh()'
#
# Generate the mesh of the dielectric_slab regression test (from its directory in
# test/data/regression/input, where mesh/mesh.jl is a link to this file):
#   julia -e 'include("mesh/mesh.jl"); generate_slab_mesh(filename="mesh/slab.msh", a=4.0,
#       b=1.0, z_min=-4.0, lengths=[3.75, 0.5, 3.75], n_z=[7, 4, 7], materials=[1, 2, 1])'

using Gmsh: gmsh

"""
    generate_slab_mesh(; filename, ...)

Generate a structured hex mesh of a periodic unit cell filled by layers stacked along z, for
the Floquet port S-parameters of planar structures at normal incidence, which the Fresnel
formulas give in closed form. Without a Floquet wave vector the fields are uniform across
the cell, so a few elements suffice laterally.

# Geometry

A unit cell of period `a` in x and `b` in y, from port 1 at z = `z_min` to port 2. Layer `i`
from port 1 has thickness `lengths[i]`, `n_z[i]` elements through it, and material
`materials[i]`: 1 for vacuum, 2 for the slab.

The defaults are the half-space of the silicon carbide example: 0.5 μm of vacuum in 8 layers,
then 0.25 μm of silicon carbide in 8 layers, with its surface at z = 0 and port 2 inside it.
The Floquet port takes the permittivity of the medium it lies in, so it terminates the
half-space exactly at any depth, and a short silicon carbide region suffices. Just below the
resonance the wavelength in silicon carbide is under 1 μm, so its layers are thin and
uniform: layers graded away from the surface would be too coarse where the wave, about 1 μm
deep there, still matters. At order 3 this mesh is within 3e-9 of the Fresnel S-parameters
from 2 to 22 μm, with about equal errors at the shortest wavelength (from the vacuum layers)
and just below the resonance (from the silicon carbide layers).

# Attributes

  - Volume 1: vacuum, Volume 2: slab
  - Boundary 1,2: x periodic pair
  - Boundary 3,4: y periodic pair
  - Boundary 5: port 1 (z = `z_min`), Boundary 6: port 2 (z = `z_min + sum(lengths)`)

The mesh has no units: `Model.L0` in the configuration file sets them (1e-6 for the silicon
carbide half-space in μm, 1e-2 for the regression test in cm).
"""
function generate_slab_mesh(;
    filename::AbstractString="mesh/sic_halfspace.msh",
    a::Real=1.0,
    b::Real=1.0,
    z_min::Real=-0.5,
    lengths::AbstractVector{<:Real}=[0.5, 0.25],
    n_z::AbstractVector{<:Integer}=[8, 8],
    materials::AbstractVector{<:Integer}=[1, 2],
    n_x::Integer=3,
    n_y::Integer=3,
    verbose::Integer=5,
    gui::Bool=false
)
    @assert length(lengths) == length(n_z) == length(materials) >= 1
    @assert all(m -> m in (1, 2), materials)

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    if "slab" in gmsh.model.list()
        gmsh.model.setCurrent("slab")
        gmsh.model.remove()
    end
    gmsh.model.add("slab")

    # Layers along z with shared interfaces.
    zs = z_min .+ [0.0; cumsum(lengths)]
    n_layers = length(lengths)
    layers =
        [gmsh.model.occ.addBox(0.0, 0.0, zs[i], a, b, zs[i + 1] - zs[i]) for i = 1:n_layers]
    if n_layers > 1
        gmsh.model.occ.fragment([(3, layers[1])], [(3, tag) for tag in layers[2:end]])
    end
    gmsh.model.occ.synchronize()
    layer_at(z) = clamp(searchsortedlast(zs, z), 1, n_layers)
    tol = 1.0e-6 * max(a, b, minimum(lengths))

    # Transfinite constraints for a structured hex mesh: each straight edge is divided
    # according to the axis it runs along, and in z according to its layer.
    for (_, tag) in gmsh.model.getEntities(1)
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(1, tag)
        n = if abs(x1 - x0) > tol
            n_x
        elseif abs(y1 - y0) > tol
            n_y
        else
            n_z[layer_at(0.5 * (z0 + z1))]
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
        push!(materials[layer_at(z)] == 2 ? slab : vacuum, tag)
    end
    sides = [Tuple{Float64, Int}[] for _ = 1:6]
    for (_, tag) in gmsh.model.getEntities(2)
        x, y, z = gmsh.model.occ.getCenterOfMass(2, tag)
        side = if abs(x) < tol
            1
        elseif abs(x - a) < tol
            2
        elseif abs(y) < tol
            3
        elseif abs(y - b) < tol
            4
        elseif abs(z - zs[1]) < tol
            5
        elseif abs(z - zs[end]) < tol
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
