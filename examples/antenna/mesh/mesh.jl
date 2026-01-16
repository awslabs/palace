# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script uses Gmsh to create a mesh for a dipole antenna enclosed within a
spherical boundary. A dipole antenna consists of two equal cylinders (the arms) separated by
a thin gap. In the thin gap, there is a flat rectangle that connects the cylinder and that
can is used as a lumped port.

The generated mesh contains five regions:
1. A 3D volume region (the space between the cylinders and the outer sphere)
2. A large outer spherical boundary (typically set to "absorbing" boundary conditions)
3. Two identical cylindrical conductors aligned on the z axis separated by a thin gap
4. A flat rectangle fills the gap between the cylinders along the xz plane

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
    # We will create this mesh with a top-down approach. We first create the 3D outer
    # sphere, the two 3D cylinders, and the 2D rectangle. Then, we fragment everything,
    # ensuring that the gap is surface is properly embedded into the end of the cylinders,
    # and that all 2D surfaces are properly embedded in the 3D space.
    #
    # Fragmenting produces:
    # - 3 3D entities (the domain and the interior of the cylinders)
    # - 1 2D entity for the rectangular port
    # - For each arm, 4 2D entities: the cylindrical area, one of the two surface caps, and
    #   the two fragments for the other cap. We get two fragments for the other cap because
    #   it is split in two by the rectangle
    #
    # After fragmenting, we identify the various pieces by looking at their expected
    # bounding boxes and we add them to the correct gmsh physical groups. Finally, we
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

    # Create geometry: outer sphere and two cylindrical arms.
    outer_boundary = kernel.addSphere(0, 0, 0, outer_boundary_radius)
    top_arm = kernel.addCylinder(0, 0, gap_size / 2, 0, 0, arm_length, arm_radius)
    bot_arm = kernel.addCylinder(0, 0, -gap_size / 2, 0, 0, -arm_length, arm_radius)

    # Create gap rectangle (port region) and rotate to XZ plane (we need to do this because
    # OpenCASCADE does not have a primitive to create a rectangle directly on the XZ
    # plane...).
    gap_rectangle = kernel.addRectangle(-arm_radius, -gap_size/2, 0, 2arm_radius, gap_size)
    kernel.rotate([(2, gap_rectangle)], 0, 0, 0, 1, 0, 0, π/2)

    # Fragment geometry to create mesh domain.
    kernel.fragment([(3, outer_boundary)], [(3, top_arm), (3, bot_arm), (2, gap_rectangle)])

    # Synchronize CAD operations with Gmsh model.
    kernel.synchronize()

    # Now that we have fragmented the geometry, the tags for the various entities might have
    # changed, so we find them again using their expected positions/sizes.

    # Helper functions to identify the various components.
    all_2d_entities = kernel.getEntities(2)
    all_3d_entities = kernel.getEntities(3)

    # zmin(x) > 0 means that x is fully above xy plane (similar with zmax(x) < 0).

    find_2D_top_arm() = filter(x -> zmin(x) > 0, all_2d_entities)
    find_2D_bot_arm() = filter(x -> zmax(x) < 0, all_2d_entities)

    # The only objects that span the entire domain have xmin = - outer_boundary_radius (this
    # is true for all the directions).
    spans_domain(x) =
        isapprox(xmin(x), -outer_boundary_radius, atol=outer_boundary_radius/100)
    find_2D_outer_sphere() = filter(spans_domain, all_2d_entities)
    find_3D_domain() = filter(spans_domain, all_3d_entities)

    # Helper functions to find 3D arm domains (interior of cylinders)
    find_3D_top_arm() = filter(x -> zmin(x) > 0 && !spans_domain(x), all_3d_entities)
    find_3D_bot_arm() = filter(x -> zmax(x) < 0 && !spans_domain(x), all_3d_entities)


    # We find the rectangle because we know its precise location and size (filling the
    # entire gap and being oriented on the XZ plane at y=0).
    function is_gap_rectangle(x)
        eps = gap_size/100
        return isapprox(xmin(x), -arm_radius, atol=eps) &&
               isapprox(xmax(x), arm_radius, atol=eps) &&
               isapprox(ymin(x), 0, atol=eps) &&
               isapprox(ymax(x), 0, atol=eps) &&
               isapprox(zmin(x), -gap_size/2, atol=eps) &&
               isapprox(zmax(x), +gap_size/2, atol=eps)
    end

    find_gap_rectangle() = filter(is_gap_rectangle, all_2d_entities)

    # Now we can identify the regions.
    outer_sphere_dimtags = find_2D_outer_sphere()
    top_arm_dimtags = find_2D_top_arm()
    bot_arm_dimtags = find_2D_bot_arm()
    domain_dimtags = find_3D_domain()
    gap_rectangle_dimtags = find_gap_rectangle()
    top_arm_domain_dimtags = find_3D_top_arm()
    bot_arm_domain_dimtags = find_3D_bot_arm()

    # Verify we found the expected number of entities.
    @assert length(top_arm_dimtags) == 4         # Cylinder has 4 surfaces (curved + top cap + 2 fragment of bottom cap)
    @assert length(bot_arm_dimtags) == 4         # Cylinder has 4 surfaces (curved + bottom cap + 2 fragment of top cap)
    @assert length(gap_rectangle_dimtags) == 1   # Single rectangular port
    @assert length(outer_sphere_dimtags) == 1    # Single outer boundary
    @assert length(top_arm_domain_dimtags) == 1  # Single 3D domain top bottom arm
    @assert length(bot_arm_domain_dimtags) == 1  # Single 3D domain top bottom arm
    @assert length(domain_dimtags) == 1          # Single 3D domain (top and bottom arms subtracted)



    # Create physical groups (these become attributes in Palace).
    top_arm_group =
        gmsh.model.addPhysicalGroup(2, extract_tag.(top_arm_dimtags), -1, "top_arm")
    bot_arm_group =
        gmsh.model.addPhysicalGroup(2, extract_tag.(bot_arm_dimtags), -1, "bot_arm")
    gap_rectangle_group =
        gmsh.model.addPhysicalGroup(2, extract_tag.(gap_rectangle_dimtags), -1, "port")
    outer_boundary_group = gmsh.model.addPhysicalGroup(
        2,
        extract_tag.(outer_sphere_dimtags),
        -1,
        "outer_boundary"
    )
    top_arm_domain_group =
        gmsh.model.addPhysicalGroup(3, extract_tag.(top_arm_domain_dimtags), -1, "top_arm_domain")
    bot_arm_domain_group =
        gmsh.model.addPhysicalGroup(3, extract_tag.(bot_arm_domain_dimtags), -1, "bot_arm_domain")
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

    # First, create a mesh size field (with id 1) that extends from the inner surfaces to
    # the volume toward the outer boundary.
    gmsh.model.mesh.field.add("Extend", 1)
    gmsh.model.mesh.field.setNumbers(
        1,
        "SurfacesList",
        extract_tag.(vcat(top_arm_dimtags, bot_arm_dimtags, gap_rectangle_dimtags))
    )
    gmsh.model.mesh.field.setNumber(1, "DistMax", outer_boundary_radius)
    gmsh.model.mesh.field.setNumber(1, "SizeMax", wavelength / n_farfield)

    # Finally, use this Extend field to determine element sizes.
    gmsh.model.mesh.field.setAsBackgroundMesh(1)

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
    println("Top arm (2D): ", top_arm_group)
    println("Bottom arm (2D): ", bot_arm_group)
    println("Rectangular port (2D): ", gap_rectangle_group)
    println("Farfield boundary (2D): ", outer_boundary_group)
    println("Top arm domain (3D): ", top_arm_domain_group)
    println("Bottom arm domain (3D): ", bot_arm_domain_group)
    println("Domain (3D): ", domain_group)
    println()

    # Optionally launch the Gmsh GUI.
    if gui
        gmsh.fltk.run()
    end

    # Clean up Gmsh resources.
    return gmsh.finalize()
end
