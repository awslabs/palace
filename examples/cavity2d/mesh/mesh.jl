# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script uses Gmsh to create a 2D mesh for a rectangular domain. The mesh is used
for 2D eigenmode, driven, and transient Palace examples.

The generated mesh contains three distinct regions:
1. A 2D surface region (the rectangular domain)
2. A bottom edge boundary (used as lumped port for driven/transient simulations)
3. The remaining edges (left, top, right — used as PEC boundaries)

## Prerequisites

This script requires the Gmsh Julia package. If you don't already have it installed, you can
install it with

```bash
julia -e 'using Pkg; Pkg.add("Gmsh")'
```

## How to run

From this directory, run:
```bash
julia -e 'include("mesh.jl"); generate_cavity2d_mesh(; filename="cavity2d.msh")'
```
=#

using Gmsh: gmsh

"""
    generate_cavity2d_mesh(;
        filename::AbstractString,
        a::Real = 1.0,
        b::Real = 0.5,
        n_elem::Integer = 10,
        verbose::Integer = 5,
        gui::Bool = false
    )

Generate a 2D rectangular mesh for the cavity2d example using Gmsh.

# Arguments

  - filename - the filename to use for the generated mesh
  - a - width of the rectangle (x direction)
  - b - height of the rectangle (y direction)
  - n_elem - approximate number of elements along the shorter side
  - verbose - flag to dictate the level of print to REPL (0-5)
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_cavity2d_mesh(;
    filename::AbstractString,
    a::Real=1.0,
    b::Real=0.5,
    n_elem::Integer=10,
    verbose::Integer=5,
    gui::Bool=false
)
    kernel = gmsh.model.occ
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    if "cavity2d" in gmsh.model.list()
        gmsh.model.setCurrent("cavity2d")
        gmsh.model.remove()
    end
    gmsh.model.add("cavity2d")

    # Create a rectangle: origin at (0,0), width a, height b
    rect = kernel.addRectangle(0.0, 0.0, 0.0, a, b)
    kernel.synchronize()

    # Get the boundary curves
    boundary = gmsh.model.getBoundary([(2, rect)], false, false, false)

    # Identify the bottom edge (y ≈ 0) vs the rest
    bottom_curves = Int[]
    other_curves = Int[]
    eps = 1.0e-6
    for (dim, tag) in boundary
        bb = gmsh.model.getBoundingBox(dim, tag)
        # bb = (xmin, ymin, zmin, xmax, ymax, zmax)
        if abs(bb[2]) < eps && abs(bb[5]) < eps
            # Both ymin and ymax are ≈ 0 → bottom edge
            push!(bottom_curves, tag)
        else
            push!(other_curves, tag)
        end
    end
    @assert length(bottom_curves) == 1 "Expected exactly one bottom edge, got $(length(bottom_curves))"
    @assert length(other_curves) == 3 "Expected three other edges, got $(length(other_curves))"

    # Create physical groups
    domain_group = gmsh.model.addPhysicalGroup(2, [rect], -1, "domain")
    bottom_group = gmsh.model.addPhysicalGroup(1, bottom_curves, -1, "bottom")
    walls_group = gmsh.model.addPhysicalGroup(1, other_curves, -1, "walls")

    # Mesh settings
    lc = b / n_elem
    gmsh.option.setNumber("Mesh.MeshSizeMin", lc / 2.0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
    gmsh.option.setNumber("Mesh.Algorithm", 6)

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(2)

    # Write in Gmsh 2.2 binary format as required by Palace
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, filename))

    println("\nFinished generating mesh. Physical group tags:")
    println("Domain: ", domain_group)
    println("Bottom edge (port): ", bottom_group)
    println("Walls (PEC): ", walls_group)
    println()

    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
