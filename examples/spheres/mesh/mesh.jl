# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script uses Gmsh to create a mesh for two conducting spheres enclosed within a
larger boundary sphere.

The generated mesh contains four distinct regions:
1. A 3D volume region (the space between spheres)
2. A large outer spherical boundary (typically set to ground potential)
3. Two inner spherical conductors (typically set as terminals)

## Prerequisites

This script requires the Gmsh Julia package. If you don't already have it installed, you can
install it with

```bash
julia -e 'using Pkg; Pkg.add("Gmsh")'
```

## How to run

From this directory, run:
```bash
julia -e 'include("mesh.jl"); generate_spheres_mesh(; filename="spheres.msh")'
```
This generates the mesh used in the example.

To visualize the mesh in Gmsh's graphical interface, add the `gui=true` parameter:
```bash
julia -e 'include("mesh.jl"); generate_spheres_mesh(; filename="spheres.msh", gui=true)'
```

The script will generate a mesh file and print the "attribute" numbers for each region.
These attributes are needed when configuring Palace simulations.
=#

using Gmsh: gmsh

# Convenience function to extract the tag (ID) from the return values of several Gmsh
# functions
extract_tag(entity) = first(entity)[2]

"""
    generate_spheres_mesh(;
        filename::AbstractString,
        radius_a::Real = 1.0,
        radius_b::Real = 2.0,
        center_d::Real = 5.0,
        radius_farfield::Real = 15.0 * center_d,
        verbose::Integer = 5,
        gui::Bool = false
    )

Generate a mesh for the two spheres example using Gmsh.

# Arguments

  - filename - the filename to use for the generated mesh
  - radius_a - radius of the first sphere
  - radius_b - radius of the second sphere
  - center_d - distance between sphere centers
  - radius_farfield - Radius of the outer boundary sphere
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
    (0-5, higher = more verbose)
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_spheres_mesh(;
    filename::AbstractString,
    radius_a::Real=1.0,
    radius_b::Real=2.0,
    center_d::Real=5.0,
    radius_farfield::Real=15.0 * center_d,
    verbose::Integer=5,
    gui::Bool=false
)
    # Boilerplate
    kernel = gmsh.model.occ
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Create a new model. The name spheres is not important. If a model was already added,
    # remove it first (this is useful when interactively the body of this function in the
    # REPL)
    if "spheres" in gmsh.model.list()
        gmsh.model.setCurrent("spheres")
        gmsh.model.remove()
    end
    gmsh.model.add("spheres")

    # n_sphere: Controls how many elements around sphere's circumference,
    #           higher values = finer mesh on spheres
    # l_farfield: Maximum element size at the outer boundary,
    #             larger values = coarser mesh at boundary
    n_sphere = 16
    l_farfield = 20.0

    # Create three spheres: two inner conductors along the x axis and one outer boundary
    # with centers along the x axis
    sphere_a = kernel.addSphere(-0.5 * center_d, 0.0, 0.0, radius_a)
    sphere_b = kernel.addSphere(0.5 * center_d, 0.0, 0.0, radius_b)
    sphere_farfield = kernel.addSphere(0.0, 0.0, 0.0, radius_farfield)

    # We want to mesh the volume between the spheres, so we start with the large outer
    # sphere and subtract the two inner spheres
    cut_dimtags, cut_dimtags_map =
        kernel.cut([(3, sphere_farfield)], [(3, sphere_a), (3, sphere_b)])

    # Verify we got exactly one 3D region as expected
    only_one_region = length(cut_dimtags) == 1
    cut_region_is_3d = first(cut_dimtags)[1] == 3
    @assert only_one_region && cut_region_is_3d

    domain = extract_tag(cut_dimtags)

    # After boolean operations, we need to find the surfaces again
    eps = 1.0e-3

    # Find sphere_a's boundary surface
    sphere_a = kernel.getEntitiesInBoundingBox(
        -0.5 * center_d - radius_a - eps,
        -radius_a - eps,
        -radius_a - eps,
        -0.5 * center_d + radius_a + eps,
        radius_a + eps,
        radius_a + eps,
        2  # 2 means we're looking for 2D surfaces
    )

    # Find sphere_b's boundary surface
    sphere_b = kernel.getEntitiesInBoundingBox(
        0.5 * center_d - radius_b - eps,
        -radius_b - eps,
        -radius_b - eps,
        0.5 * center_d + radius_b + eps,
        radius_b + eps,
        radius_b + eps,
        2
    )
    @assert length(sphere_a) == 1 && length(sphere_b) == 1

    sphere_a = extract_tag(sphere_a)
    sphere_b = extract_tag(sphere_b)

    # Find the outer boundary by getting all 2D entities except the inner spheres
    sphere_farfield =
        filter(x -> x != (2, sphere_a) && x != (2, sphere_b), kernel.getEntities(2))
    @assert length(sphere_farfield) == 1
    sphere_farfield = extract_tag(sphere_farfield)

    # Commit all geometric operations from the CAD representation to the Gmsh model
    kernel.synchronize()

    # Create physical groups (these become attributes in Palace)
    # The provided names are primarily for human consumption, e.g., in the Gmsh GUI
    # The -1 means "assign an attribute number automatically"
    domain_group = gmsh.model.addPhysicalGroup(3, [domain], -1, "domain")
    farfield_group = gmsh.model.addPhysicalGroup(2, [sphere_farfield], -1, "farfield")
    sphere_a_group = gmsh.model.addPhysicalGroup(2, [sphere_a], -1, "sphere_a")
    sphere_b_group = gmsh.model.addPhysicalGroup(2, [sphere_b], -1, "sphere_b")

    # Set minimum nodes per curve for better curved element quality
    gmsh.option.setNumber("Mesh.MinimumCurveNodes", 2)
    gmsh.option.setNumber("Mesh.MinimumCircleNodes", 0)

    # Set smallest and largest element size
    gmsh.option.setNumber("Mesh.MeshSizeMin", 2.0 * pi * radius_a / n_sphere / 2.0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
    # Set minimum number of elements per 2π radians of curvature
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", n_sphere)
    # Don't extend mesh size constraints from boundaries into the volume
    # This option is typically activated when working with mesh size fields
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    # Create mesh size fields to manually control size and distribution of elements
    # throughout the mesh

    # First, create a mesh size field (with id 1) that extends from the inner surfaces to
    # the volume toward the outer boundary
    gmsh.model.mesh.field.add("Extend", 1)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", [sphere_a, sphere_b])
    gmsh.model.mesh.field.setNumber(1, "DistMax", radius_farfield)
    gmsh.model.mesh.field.setNumber(1, "SizeMax", l_farfield)

    # Finally, use this Extend field to determine element sizes
    gmsh.model.mesh.field.setAsBackgroundMesh(1)

    # Choose meshing algorithm. Typically, we would choose HXT, 10, because it
    # is parallel and high-performance, but it is not reproducible, so best to
    # stick with something more stable for this example.
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    gmsh.model.mesh.generate(3) # 3 means generate a 3D volume mesh
    gmsh.model.mesh.setOrder(3) # 3 means cubically curved elements

    # Set mesh format version as required by Palace
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, filename))

    println("\nFinished generating mesh. Physical group tags:")
    println("Domain: ", domain_group)
    println("Farfield boundary: ", farfield_group)
    println("Sphere A: ", sphere_a_group)
    println("Sphere B: ", sphere_b_group)
    println()

    # Optionally launch the Gmsh GUI
    if gui
        gmsh.fltk.run()
    end

    # Clean up Gmsh resources
    return gmsh.finalize()
end
