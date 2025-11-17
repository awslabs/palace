# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script uses Gmsh to create a mesh for a dipole antenna enclosed within a
spherical boundary. A dipole antenna consists of two equal cylinders (the arms) separated by
a thin gap. In the thin gap, there is a flat rectangle that connects the cylinder and that
can is used as a lumped port.

The generated mesh contains two regions:
1. A 3D volume region (the space inside the sphere)
2. A large outer spherical boundary (typically set to "absorbing" boundary conditions)

## Prerequisites

This script requires the Gmsh Julia package. If you don't already have it installed, you can
install it with

```bash
julia -e 'using Pkg; Pkg.add("Gmsh")'
```

## How to run

From this directory, run:
```bash
julia -e 'include("mesh.jl"); generate_antenna_mesh(; filename="antenna.msh")'
```
To visualize the mesh in Gmsh's graphical interface, add the `gui=true` parameter:
```bash
julia -e 'include("mesh.jl"); generate_antenna_mesh(; filename="antenna.msh", gui=true)'
```

The script will generate a mesh file and print the "attribute" numbers for each region.
These attributes are needed when configuring Palace simulations.
=#

using Gmsh: gmsh

"""
    extract_tag(object)

Extract the Gmsh tag in `object`.

If `object` contains only one tag, return it as an integer, otherwise, preserve its
container.

Most gmsh functions return list of tuples like `[(2, 5), (2, 8), (2, 10), ...]`, where the
first number is dimensionality and the second is the integer tag associated to that object.

#### Example

```jldoctest
julia> extract_tag((3, 6))
6

julia> entities = [(2, 5), (2, 8), (2, 10)];

julia> extract_tag.(entities)
[5, 8, 10]
```
"""
extract_tag(object) = extract_tag(only(object))
extract_tag(object::Tuple) = last(object)
extract_tag(object::Integer) = error("You passed an integer tag directly to `extract_tag`")

# Convenience functions to extract extrema of the bounding box of an entity
xmin(x::Tuple) = gmsh.model.occ.get_bounding_box(x...)[1]
ymin(x::Tuple) = gmsh.model.occ.get_bounding_box(x...)[2]
zmin(x::Tuple) = gmsh.model.occ.get_bounding_box(x...)[3]
xmax(x::Tuple) = gmsh.model.occ.get_bounding_box(x...)[4]
ymax(x::Tuple) = gmsh.model.occ.get_bounding_box(x...)[5]
zmax(x::Tuple) = gmsh.model.occ.get_bounding_box(x...)[6]

"""
    generate_antenna_mesh(;
                          filename::AbstractString,
                          wavelength::Real=4.0,
                          arm_length::Real=wavelength/4,
                          arm_radius::Real=arm_length/20,
                          gap_size::Real=arm_length/100,
                          outer_boundary_radius::Real=1.5wavelength,
                          verbose::Integer=5,
                          gui::Bool=false
                         )

Generate a mesh for a dipole antenna using Gmsh.

# Arguments

  - filename - output mesh filename
  - wavelength - wavelength of the resulting electromagnetic wave
  - arm_length - length of each antenna arm
  - arm_radius - radius of the cylindrical antenna arms
  - gap_size - size of the gap between the two arms (port region)
  - outer_boundary_radius - radius of the outer spherical boundary
  - verbose -  gmsh verbosity level (0-5, higher = more verbose)
  - gui - whether to launch the Gmsh GUI after mesh generation
"""
function generate_antenna_mesh(;
    filename::AbstractString,
    wavelength::Real=4.0,
    arm_length::Real=wavelength/4,
    arm_radius::Real=arm_length/20,
    gap_size::Real=arm_length/100,
    outer_boundary_radius::Real=1.5wavelength,
    verbose::Integer=5,
    gui::Bool=false
)
    # We will create this mesh with a simple approach. We create the 3D
    # sphere only, which produces:
    # - 1 3D entity (the domain)
    # - 1 2D entity (the outer boundary)
    #
    # After creating, we add them to the correct gmsh physical groups. Finally, we
    # control mesh size with a mesh size field and generate the mesh.

    # Boilerplate
    gmsh.initialize()
    kernel = gmsh.model.occ
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Create a new model. The name dipole is not important. If a model was already added,
    # remove it first (this is useful when interactively evaluating the body of this
    # function in the REPL).
    if "dipole" in gmsh.model.list()
        gmsh.model.setCurrent("dipole")
        gmsh.model.remove()
    end
    gmsh.model.add("dipole")

    # Mesh refinement parameter: controls elements around cylinder circumference.
    # Higher number = higher resolution.
    n_circle = 12
    # How many elements per wavelength on the outer sphere.
    # Higher number = higher resolution.
    n_farfield = 3

    # Create geometry
    outer_boundary = kernel.addSphere(0, 0, 0, outer_boundary_radius)

    # Synchronize CAD operations with Gmsh model.
    kernel.synchronize()

    # Helper functions to identify the various components.
    all_2d_entities = kernel.getEntities(2)
    all_3d_entities = kernel.getEntities(3)

    # For a simple sphere, there should be exactly one 2D and one 3D entity.
    outer_sphere_dimtags = all_2d_entities
    domain_dimtags = all_3d_entities

    # Verify we found the expected number of entities.
    @assert length(outer_sphere_dimtags) == 1  # Single outer boundary
    @assert length(domain_dimtags) == 1        # Single 3D domain

    # Create physical groups (these become attributes in Palace).
    outer_boundary_group = gmsh.model.addPhysicalGroup(
        2,
        extract_tag.(outer_sphere_dimtags),
        -1,
        "outer_boundary"
    )
    domain_group =
        gmsh.model.addPhysicalGroup(3, extract_tag.(domain_dimtags), -1, "domain")

    # Set mesh size parameters.
    gmsh.option.setNumber("Mesh.MeshSizeMin", 2.0 * pi * arm_radius / n_circle / 2.0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", wavelength / n_farfield)
    # Set minimum number of elements per 2π radians of curvature.
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", n_circle)
    # Don't extend mesh size constraints from boundaries into the volume.
    # This option is typically activated when working with mesh size fields.
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    # Finally, we control mesh size using a mesh size field.

    # Create a simple distance field from the center for mesh sizing.
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "PointsList", [])  # No specific points

    # Use a Box field to control mesh size - finer at center, coarser toward boundary
    gmsh.model.mesh.field.add("Box", 2)
    gmsh.model.mesh.field.setNumber(2, "VIn", wavelength / n_farfield / 2)  # Fine mesh at center
    gmsh.model.mesh.field.setNumber(2, "VOut", wavelength / n_farfield)     # Coarser toward boundary
    gmsh.model.mesh.field.setNumber(2, "XMin", -outer_boundary_radius/4)
    gmsh.model.mesh.field.setNumber(2, "XMax", outer_boundary_radius/4)
    gmsh.model.mesh.field.setNumber(2, "YMin", -outer_boundary_radius/4)
    gmsh.model.mesh.field.setNumber(2, "YMax", outer_boundary_radius/4)
    gmsh.model.mesh.field.setNumber(2, "ZMin", -outer_boundary_radius/4)
    gmsh.model.mesh.field.setNumber(2, "ZMax", outer_boundary_radius/4)

    # Use this Box field to determine element sizes.
    gmsh.model.mesh.field.setAsBackgroundMesh(2)

    # Set 2D/3D meshing algorithm. Chosen to be deterministic, not necessarily the
    # best.
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)

    # Generate 3D volume mesh and set to 3rd order elements.
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(3)

    # Set output format for Palace compatibility.
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, filename))

    println("\nFinished generating mesh. Physical group tags:")
    println("Farfield boundary (2D): ", outer_boundary_group)
    println("Domain (3D): ", domain_group)
    println()

    # Optionally launch the Gmsh GUI.
    if gui
        gmsh.fltk.run()
    end

    # Clean up Gmsh resources.
    return gmsh.finalize()
end
