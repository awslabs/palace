# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Gmsh: gmsh

"""
    generate_cylindrical_cavity_mesh(;
        refinement=0,
        order=1,
        mesh_type=0,
        radius=2.74,
        aspect_ratio=1.0,
        filename,
        verbose=1,
    )

Generate a mesh for the cylindrical cavity resonator example using Gmsh

# Arguments

  - refinement - measure of how many elements to include, 0 is least
  - order - the polynomial order of the approximation, minimum 1
  - mesh_type - 0 = tetrahedral mesh, 1 = prism mesh, 2 = hexahedral mesh
  - radius - the radius of the cavity resonator
  - aspect_ratio - the ratio of the DIAMETER of the cavity to the height
  - filename - the filename to use for the generated mesh
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
"""
function generate_cylindrical_cavity_mesh(;
    refinement::Integer=0,
    order::Integer=1,
    mesh_type::Integer=0,
    radius::Real=2.74,
    aspect_ratio::Real=1.0,
    filename::AbstractString,
    verbose::Integer=1
)
    @assert refinement >= 0
    @assert order > 0

    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "cavity" in gmsh.model.list()
        gmsh.model.setCurrent("cavity")
        gmsh.model.remove()
    end
    gmsh.model.add("cavity")

    # Geometry parameters (in cm)
    height = aspect_ratio * 2 * radius  # Cylinder height

    # Mesh parameters
    n_height = 2                       # Two elements in height
    n_circum = mesh_type == 2 ? 4 : 6  # Four or six elements on circumference

    # Geometry
    if (mesh_type == 2)
        base_square = kernel.addRectangle(
            -0.4 * radius,
            -0.4 * radius,
            0.0,
            0.8 * radius,
            0.8 * radius
        )
        kernel.rotate((2, base_square), 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, pi / 4.0)
        base_circle = kernel.addDisk(0.0, 0.0, 0.0, radius, radius)
        base_circle, _ = kernel.cut((2, base_circle), (2, base_square), -1, true, false)
        @assert length(base_circle) == 1 && first(base_circle)[1] == 2
        base_circle = first(base_circle)[2]
        cylinder_dimtags = kernel.extrude(
            [(2, base_square), (2, base_circle)],
            0.0,
            0.0,
            height,
            [n_height],
            [1.0],
            true
        )
        cylinder = filter(x -> x[1] == 3, cylinder_dimtags)
        @assert length(cylinder) == 2
        cylinder = [x[2] for x in cylinder]

        kernel.synchronize()

        gmsh.model.mesh.setRecombine(2, base_square)
        gmsh.model.mesh.setRecombine(2, base_circle)

        # Add physical groups
        cylinder_group = gmsh.model.addPhysicalGroup(3, cylinder, -1, "cylinder")
        _, square_boundaries = gmsh.model.getAdjacencies(3, cylinder[1])
        _, circle_boundaries = gmsh.model.getAdjacencies(3, cylinder[2])
        boundaries = copy(circle_boundaries)
        for boundary in square_boundaries
            idx = first(indexin([boundary], boundaries))
            if !isnothing(idx)
                deleteat!(boundaries, idx)
            else
                push!(boundaries, boundary)
            end
        end
        boundary_group = gmsh.model.addPhysicalGroup(2, boundaries, -1, "boundaries")
    else
        base_circle = kernel.addDisk(0.0, 0.0, 0.0, radius, radius)
        cylinder_dimtags = kernel.extrude(
            [(2, base_circle)],
            0.0,
            0.0,
            height,
            [n_height],
            [1.0],
            mesh_type > 0
        )
        cylinder = filter(x -> x[1] == 3, cylinder_dimtags)
        @assert length(cylinder) == 1 && first(cylinder)[1] == 3
        cylinder = first(cylinder)[2]

        kernel.synchronize()

        # Add physical groups
        cylinder_group = gmsh.model.addPhysicalGroup(3, [cylinder], -1, "cylinder")

        _, boundaries = gmsh.model.getAdjacencies(3, cylinder)
        boundary_group = gmsh.model.addPhysicalGroup(2, boundaries, -1, "boundaries")
    end

    # Generate mesh
    gmsh.option.setNumber("Mesh.MinimumCurveNodes", 2)
    gmsh.option.setNumber("Mesh.MinimumCircleNodes", 0)

    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", n_circum)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)

    gmsh.option.setNumber("Mesh.Algorithm", mesh_type == 2 ? 8 : 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)

    gmsh.model.mesh.generate(3) # Dimension of the mesh
    for i = 0:(refinement - 1)
        gmsh.model.mesh.refine()
    end
    gmsh.model.mesh.setOrder(order) # Polynomial order of the mesh

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(joinpath(@__DIR__, filename))

    # Print some information
    if verbose > 0
        println("\nFinished generating mesh. Physical group tags:")
        println("Cylinder: ", cylinder_group)
        println("Boundaries: ", boundary_group)
        println()
    end

    # Optionally launch GUI
    if "gui" in lowercase.(ARGS)
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
