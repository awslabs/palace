# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script uses Gmsh to create a 2D mesh for an annular domain (two concentric
circles). The mesh is used for 2D electrostatic and magnetostatic Palace examples.

The generated mesh contains three distinct regions:
1. A 2D surface region (the annular domain between inner and outer circles)
2. An inner circular boundary (terminal or surface current)
3. An outer circular boundary (ground or PEC)

## Prerequisites

This script requires the Gmsh Julia package. If you don't already have it installed, you can
install it with

```bash
julia -e 'using Pkg; Pkg.add("Gmsh")'
```

## How to run

From this directory, run:
```bash
julia -e 'include("mesh.jl"); generate_disc2d_mesh(; filename="disc2d.msh")'
```
=#

using Gmsh: gmsh

"""
    generate_disc2d_mesh(;
        filename::AbstractString,
        radius_inner::Real = 1.0,
        radius_outer::Real = 4.0,
        n_circle::Integer = 24,
        verbose::Integer = 5,
        gui::Bool = false
    )

Generate a 2D annular mesh for the disc2d example using Gmsh.

# Arguments

  - filename - the filename to use for the generated mesh
  - radius_inner - radius of the inner circle
  - radius_outer - radius of the outer circle
  - n_circle - minimum number of elements per 2π radians of curvature
  - verbose - flag to dictate the level of print to REPL (0-5)
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_disc2d_mesh(;
    filename::AbstractString,
    radius_inner::Real=1.0,
    radius_outer::Real=4.0,
    n_circle::Integer=24,
    verbose::Integer=5,
    gui::Bool=false
)
    kernel = gmsh.model.occ
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    if "disc2d" in gmsh.model.list()
        gmsh.model.setCurrent("disc2d")
        gmsh.model.remove()
    end
    gmsh.model.add("disc2d")

    # Create two concentric disks and subtract the inner from the outer
    disk_outer = kernel.addDisk(0.0, 0.0, 0.0, radius_outer, radius_outer)
    disk_inner = kernel.addDisk(0.0, 0.0, 0.0, radius_inner, radius_inner)
    cut_dimtags, _ = kernel.cut([(2, disk_outer)], [(2, disk_inner)])
    @assert length(cut_dimtags) == 1 && cut_dimtags[1][1] == 2

    domain = cut_dimtags[1][2]

    # Get boundary curves
    kernel.synchronize()
    boundary = gmsh.model.getBoundary([(2, domain)], false, false, false)
    @assert length(boundary) == 2

    # Identify inner vs. outer circle by checking a point on the curve
    inner_curve = -1
    outer_curve = -1
    for (dim, tag) in boundary
        bb = gmsh.model.getBoundingBox(dim, tag)
        # bb = (xmin, ymin, zmin, xmax, ymax, zmax)
        extent = max(bb[4] - bb[1], bb[5] - bb[2])
        if extent < 2.0 * (radius_inner + radius_outer) / 2.0
            inner_curve = tag
        else
            outer_curve = tag
        end
    end
    @assert inner_curve > 0 && outer_curve > 0

    # Create physical groups (these become attributes in Palace)
    domain_group = gmsh.model.addPhysicalGroup(2, [domain], -1, "domain")
    inner_group = gmsh.model.addPhysicalGroup(1, [inner_curve], -1, "inner")
    outer_group = gmsh.model.addPhysicalGroup(1, [outer_curve], -1, "outer")

    # Mesh settings
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", n_circle)
    gmsh.option.setNumber("Mesh.MeshSizeMin", 2.0 * pi * radius_inner / n_circle / 2.0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 2.0 * pi * radius_outer / n_circle)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(2)

    # Write in Gmsh 2.2 binary format as required by Palace
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, filename))

    println("\nFinished generating mesh. Physical group tags:")
    println("Domain: ", domain_group)
    println("Inner boundary: ", inner_group)
    println("Outer boundary: ", outer_group)
    println()

    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
