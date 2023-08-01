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
    verbose::Integer=1,
    gui::Bool=false,
    internal_boundary::Bool=false
)
    @assert refinement >= 0
    @assert order > 0

    kernel = gmsh.model.occ

    gmsh.initialize()
    # gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "cavity" in gmsh.model.list()
        gmsh.model.setCurrent("cavity")
        gmsh.model.remove()
    end
    gmsh.model.add("cavity")

    # Geometry parameters (in cm)
    height = aspect_ratio * 2 * radius  # Cylinder height

    # Mesh parameters
    n_height = 2 * 2^refinement # Minimum two elements in vertical
    n_circum = 8 * 2^refinement # Minimum four elements on round

    # Geometry
    base_circle = kernel.addDisk(0.0, 0.0, 0.0, radius, radius)
    if mesh_type > 0
        cylinder_dimtags =
            kernel.extrude([(2, base_circle)], 0.0, 0.0, height, [n_height], [1.0], true)
    else
        cylinder_dimtags = kernel.extrude([(2, base_circle)], 0.0, 0.0, height)
    end
    cylinder = filter(x -> x[1] == 3, cylinder_dimtags)
    @assert length(cylinder) == 1 && first(cylinder)[1] == 3
    cylinder = first(cylinder)[2]

    kernel.synchronize()

    # Add physical groups
    cylinder_group = gmsh.model.addPhysicalGroup(3, [cylinder], -1, "cylinder")

    _, boundaries = gmsh.model.getAdjacencies(3, cylinder)
    boundary_group = gmsh.model.addPhysicalGroup(2, boundaries, -1, "boundaries")

    # Generate mesh
    gmsh.option.setNumber("Mesh.MinimumCurveNodes", 4)
    gmsh.option.setNumber("Mesh.MinimumCircleNodes", 4)

    gmsh.option.setNumber("Mesh.MeshSizeMin", 2π * radius / n_circum)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 2π * radius / n_circum)

    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", n_circum)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
    if mesh_type > 1
        gmsh.model.mesh.setRecombine(2, base_circle)
    end

    # gmsh.option.setNumber("Mesh.Algorithm", 6)
    # gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    gmsh.option.setNumber("Mesh.Smoothing", 100)

    gmsh.model.mesh.generate(3) # Dimension of the mesh
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
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
