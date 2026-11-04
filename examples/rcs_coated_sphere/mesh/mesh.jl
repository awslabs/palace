# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script uses Gmsh to create the mesh of a plane-wave scattering problem with the
total-field/scattered-field (TF/SF) formulation: a coated sphere (a PEC core with a
dielectric or magnetic coating) at the center of three nested cubes,

  - the total-field cube, whose surface S is the Huygens surface of the incident plane
    wave (the field is the total field inside S and the scattered field outside);
  - the far-field cube, an interior surface outside S on which the far field of the
    scattered field is computed;
  - the outer box, the truncation of the domain (absorbing boundary).

The cubes are meshed with tetrahedra, with second-order (curved) elements to resolve the
spherical surfaces. With `empty=true` there is no sphere, which measures the leakage of the
TF/SF source (the scattered field of an empty box vanishes).

The generated mesh contains the following attributes:

Volumes:
  - 1: coating (r_core < r < r_coat)
  - 2: vacuum inside the total-field cube (half-side l_tf)
  - 3: vacuum outside the total-field cube (scattered field)

Boundaries:
  - 10: surface of the PEC core (r = r_core)
  - 11: far-field cube (half-side l_ff), interior surface
  - 12: outer boundary (half-side l_box)
  - 13: Huygens surface S (half-side l_tf), interior surface

## Prerequisites

This script requires the Gmsh Julia package. If you don't already have it installed, you can
install it with

```bash
julia -e 'using Pkg; Pkg.add("Gmsh")'
```

## How to run

From this directory, run:
```bash
julia -e 'include("mesh.jl"); generate_rcs_coated_sphere_mesh(; filename="rcs_coated_sphere.msh"); generate_rcs_coated_sphere_mesh(; filename="rcs_empty.msh", empty=true)'
```
To visualize the mesh in Gmsh's graphical interface, add the `gui=true` parameter.
=#

using Gmsh: gmsh

"""
    generate_rcs_coated_sphere_mesh(;
        filename::AbstractString,
        r_core::Real   = 0.05,
        r_coat::Real   = 0.075,
        l_tf::Real     = 0.1,
        l_ff::Real     = 0.18,
        l_box::Real    = 0.3,
        h::Real        = 0.05,
        h_coat::Real   = 0.0125,
        empty::Bool    = false,
        order::Integer = 2,
        verbose::Integer = 5,
        gui::Bool      = false
    )

Generate a mesh of a coated sphere inside the total-field, far-field and outer cubes.

# Arguments

  - filename - the filename to use for the generated mesh
  - r_core - radius of the PEC core (a hole in the mesh), in m
  - r_coat - outer radius of the coating, in m
  - l_tf - half-side of the total-field cube (Huygens surface S), in m
  - l_ff - half-side of the far-field cube, in m
  - l_box - half-side of the outer box, in m
  - h - element size in vacuum, in m
  - h_coat - element size in the coating, in m
  - empty - no sphere, only the cubes (to measure the leakage of the TF/SF source)
  - order - polynomial order of the geometry
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI after mesh generation
"""
function generate_rcs_coated_sphere_mesh(;
    filename::AbstractString,
    r_core::Real=0.05,
    r_coat::Real=0.075,
    l_tf::Real=0.1,
    l_ff::Real=0.18,
    l_box::Real=0.3,
    h::Real=0.05,
    h_coat::Real=0.0125,
    empty::Bool=false,
    order::Integer=2,
    verbose::Integer=5,
    gui::Bool=false
)
    @assert empty || 0 < r_core < r_coat < l_tf
    # The far-field cube must be outside of the elements which touch S, where the field is
    # still the total field.
    @assert l_tf + h < l_ff < l_box

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    if "rcs_coated_sphere" in gmsh.model.list()
        gmsh.model.setCurrent("rcs_coated_sphere")
        gmsh.model.remove()
    end
    gmsh.model.add("rcs_coated_sphere")
    kernel = gmsh.model.occ

    cube(a) = kernel.addBox(-a, -a, -a, 2a, 2a, 2a)
    tools = [cube(l_tf), cube(l_ff), cube(l_box)]
    if !empty
        # Coating with the PEC core as a hole.
        core = kernel.addSphere(0, 0, 0, r_core)
        ball = kernel.addSphere(0, 0, 0, r_coat)
        coat, _ = kernel.cut([(3, ball)], [(3, core)])
        pushfirst!(tools, last(only(coat)))
    end
    kernel.fragment([(3, t) for t in tools], Tuple{Int32, Int32}[])
    kernel.synchronize()

    # Classify the volumes by the half-extent of their bounding box. The fragment also
    # leaves a volume filling the core, which is removed so that the PEC core is a hole.
    halfsize(dim, tag) = maximum(diff(reshape(collect(gmsh.model.getBoundingBox(dim, tag)), 3, 2), dims=2)) / 2
    core_vols, coat_vols, tf_vols, sf_vols = Int[], Int[], Int[], Int[]
    for (_, v) in gmsh.model.getEntities(3)
        a = halfsize(3, v)
        if !empty && a <= 1.001 * r_core
            push!(core_vols, v)
        elseif !empty && a <= 1.001 * r_coat
            push!(coat_vols, v)
        elseif a <= 1.001 * l_tf
            push!(tf_vols, v)
        else
            push!(sf_vols, v)
        end
    end
    @assert length(core_vols) == (empty ? 0 : 1)
    @assert length(coat_vols) == (empty ? 0 : 1)
    if !empty
        gmsh.model.removeEntities([(3, v) for v in core_vols])
        gmsh.model.addPhysicalGroup(3, coat_vols, 1, "coating")
    end
    gmsh.model.addPhysicalGroup(3, tf_vols, 2, "total_field")
    gmsh.model.addPhysicalGroup(3, sf_vols, 3, "scattered_field")

    # Classify the surfaces: the sphere of the core by its extent, and the faces of the
    # cubes by the plane they lie on.
    tol = 1.0e-6
    core_srf, ff_srf, box_srf, tf_srf = Int[], Int[], Int[], Int[]
    for (_, s) in gmsh.model.getEntities(2)
        bb = gmsh.model.getBoundingBox(2, s)
        lo, hi = bb[1:3], bb[4:6]
        ext = maximum(abs.(bb))
        flat = minimum(hi .- lo) < tol
        if !flat && !empty && abs(ext - r_core) < 1.0e-4
            push!(core_srf, s)
        elseif flat && abs(ext - l_ff) < tol
            push!(ff_srf, s)
        elseif flat && abs(ext - l_box) < tol
            push!(box_srf, s)
        elseif flat && abs(ext - l_tf) < tol
            push!(tf_srf, s)
        end
    end
    @assert length(core_srf) == (empty ? 0 : 1)
    @assert length(ff_srf) == 6 && length(box_srf) == 6 && length(tf_srf) == 6
    if !empty
        gmsh.model.addPhysicalGroup(2, core_srf, 10, "pec_core")
    end
    gmsh.model.addPhysicalGroup(2, ff_srf, 11, "farfield")
    gmsh.model.addPhysicalGroup(2, box_srf, 12, "outer")
    gmsh.model.addPhysicalGroup(2, tf_srf, 13, "huygens")

    # Element size: h in vacuum, h_coat in the coating (and a smooth transition around it).
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 12)
    if !empty
        field = gmsh.model.mesh.field.add("Ball")
        gmsh.model.mesh.field.setNumber(field, "Radius", 1.02 * r_coat)
        gmsh.model.mesh.field.setNumber(field, "VIn", h_coat)
        gmsh.model.mesh.field.setNumber(field, "VOut", h)
        gmsh.model.mesh.field.setAsBackgroundMesh(field)
    end
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)
    if order > 1
        # Untangle the curved elements next to the spheres.
        gmsh.model.mesh.optimize("HighOrder")
    end
    _, tags, _ = gmsh.model.mesh.getElements(3)
    quality = gmsh.model.mesh.getElementQualities(collect(Iterators.flatten(tags)), "minSJ")
    @assert minimum(quality) > 0 "invalid curved elements (reduce h_coat)"

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, filename))

    println("\nFinished generating mesh. Physical group tags:")
    if !empty
        println("Coating (3D): 1")
        println("PEC core (2D): 10")
    end
    println("Total-field region (3D): 2")
    println("Scattered-field region (3D): 3")
    println("Far-field surface (2D): 11")
    println("Outer boundary (2D): 12")
    println("Huygens surface S (2D): 13")
    println("Tetrahedra: ", sum(length, tags))
    println()

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
