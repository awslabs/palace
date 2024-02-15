# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Gmsh: gmsh

"""
    generate_cylindrical_cavity_mesh(;
        filename::AbstractString,
        refinement::Integer  = 0,
        order::Integer       = 1,
        mesh_type::Integer   = 0,
        radius::Real         = 2.74,
        aspect_ratio::Real   = 1.0,
        symmetry_plane::Real = true,
        verbose::Integer     = 5,
        gui::Bool            = false
    )

Generate a mesh for the cylindrical cavity resonator example using Gmsh

# Arguments

  - filename - the filename to use for the generated mesh
  - refinement - measure of how many elements to include, 0 is least
  - order - the polynomial order of the approximation, minimum 1
  - mesh_type - 0 = tetrahedral mesh, 1 = prism mesh, 2 = hexahedral mesh
  - radius - the radius of the cavity resonator
  - aspect_ratio - the ratio of the diameter (not radius!) of the cavity to the height
  - symmetry_plane - whether to cut the cylinder in half and use a symmetry plane in the model
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_cylindrical_cavity_mesh(;
    filename::AbstractString,
    refinement::Integer  = 0,
    order::Integer       = 1,
    mesh_type::Integer   = 0,
    radius::Real         = 2.74,
    aspect_ratio::Real   = 1.0,
    symmetry_plane::Real = true,
    verbose::Integer     = 5,
    gui::Bool            = false
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
    n_circum = mesh_type == 2 ? 2 : 6  # Four or six elements on circumference

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
        base = [base_circle, base_square]
    else
        base_circle = kernel.addDisk(0.0, 0.0, 0.0, radius, radius)
        base = [base_circle]
    end

    if symmetry_plane
        cut_tool =
            kernel.addBox(-radius, 0.0, -0.1 * height, 2.0 * radius, radius, 0.1 * height)
        base_dimtags, _ = kernel.intersect([(2, x) for x in base], [(3, cut_tool)])
        base = last.(filter(x -> x[1] == 2, base_dimtags))
    end

    cylinder_dimtags = kernel.extrude(
        [(2, x) for x in base],
        0.0,
        0.0,
        height,
        [n_height],
        [1.0],
        mesh_type > 0
    )
    cylinder = last.(filter(x -> x[1] == 3, cylinder_dimtags))

    kernel.synchronize()

    # Add physical groups
    cylinder_group = gmsh.model.addPhysicalGroup(3, cylinder, -1, "cylinder")

    symmetry =
        last.(
            gmsh.model.getEntitiesInBoundingBox(
                -1.1 * radius,
                -0.1 * radius,
                -0.1 * height,
                1.1 * radius,
                0.1 * radius,
                1.1 * height,
                2
            )
        )
    boundaries = []
    for domain in cylinder
        _, domain_boundaries = gmsh.model.getAdjacencies(3, domain)
        for boundary in domain_boundaries
            if boundary in symmetry
                continue
            end
            idx = first(indexin([boundary], boundaries))
            if isnothing(idx)
                push!(boundaries, boundary)
            else
                deleteat!(boundaries, idx)
            end
        end
    end
    boundary_group = gmsh.model.addPhysicalGroup(2, boundaries, -1, "boundaries")
    symmetry_group = gmsh.model.addPhysicalGroup(2, symmetry, -1, "symmetry")

    # Generate mesh
    gmsh.option.setNumber("Mesh.MinimumCurveNodes", 2)
    gmsh.option.setNumber("Mesh.MinimumCircleNodes", 0)

    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", n_circum)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)

    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)

    if (mesh_type == 2)
        base_boundaries =
            last.(
                gmsh.model.getEntitiesInBoundingBox(
                    -1.1 * radius,
                    -1.1 * radius,
                    -0.1 * height,
                    1.1 * radius,
                    1.1 * radius,
                    0.1 * height,
                    2
                )
            )
        for boundary in base_boundaries
            gmsh.model.mesh.setRecombine(2, boundary)
        end
    end

    gmsh.model.mesh.generate(3)
    for i = 0:(refinement - 1)
        gmsh.model.mesh.refine()
    end
    gmsh.model.mesh.setOrder(order)

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(joinpath(@__DIR__, filename))

    # Print some information
    if verbose > 0
        println("\nFinished generating mesh. Physical group tags:")
        println("Cylinder: ", cylinder_group)
        println("Boundaries: ", boundary_group)
        if length(symmetry) > 0
            println("Symmetry boundaries: ", symmetry_group)
        end
    end

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
