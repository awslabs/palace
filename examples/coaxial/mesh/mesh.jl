# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Gmsh: gmsh

"""
    generate_coaxial_mesh(;
        filename::AbstractString,
        refinement::Integer     = 0,
        order::Integer          = 1,
        inner_diameter_mm::Real = 1.6383,
        outer_diameter_mm::Real = 5.461,
        length_mm::Real         = 40.0,
        verbose::Integer        = 5,
        gui::Bool               = false
    )

Generate a mesh for the coaxial cable example using Gmsh

# Arguments

  - filename - the filename to use for the generated mesh
  - refinement - measure of how many elements to include, 0 is least
  - order - the polynomial order of the approximation, minimum 1
  - inner_diameter_mm - the inner diameter of the cable, in millimeters
  - outer_diameter_mm - the outer diameter of the cable, in millimeters
  - length_mm - the length of the cable, in millimeters
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_coaxial_mesh(;
    filename::AbstractString,
    refinement::Integer     = 0,
    order::Integer          = 1,
    inner_diameter_mm::Real = 1.6383,
    outer_diameter_mm::Real = 5.461,
    length_mm::Real         = 40.0,
    verbose::Integer        = 5,
    gui::Bool               = false
)
    @assert outer_diameter_mm > inner_diameter_mm > 0
    @assert length_mm > 0
    @assert refinement >= 0
    @assert order > 0

    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "coaxial" in gmsh.model.list()
        gmsh.model.setCurrent("coaxial")
        gmsh.model.remove()
    end
    gmsh.model.add("coaxial")

    # Geometry parameters (in mm)
    ri = inner_diameter_mm / 2
    ro = outer_diameter_mm / 2

    # Mesh parameters
    n_circum = 4 * 2^refinement # min 4 elements in round
    n_length = 4 * 2^refinement # min 4 elements on length

    # Geometry
    p0 = kernel.addPoint(ri, 0.0, 0.0)
    p1 = kernel.addPoint(ro, 0.0, 0.0)

    l0 = kernel.addLine(p0, p1)

    base_face_0 = kernel.revolve(
        (1, l0),
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        pi,
        [n_circum ÷ 2],
        [1.0],
        true
    )
    base_face_1 = kernel.revolve(
        (1, l0),
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        -pi,
        [n_circum ÷ 2],
        [1.0],
        true
    )
    filter!(x -> x[1] == 2, base_face_0)
    filter!(x -> x[1] == 2, base_face_1)
    @assert length(base_face_0) == 1 && length(base_face_1) == 1

    cylinder_0 = kernel.extrude(base_face_0, 0.0, 0.0, length_mm, [n_length], [1.0], true)
    cylinder_1 = kernel.extrude(base_face_1, 0.0, 0.0, length_mm, [n_length], [1.0], true)
    base_face_0 = last(first(base_face_0))
    base_face_1 = last(first(base_face_1))
    far_face_0 = last(cylinder_0[1])
    cylinder_0 = last(cylinder_0[2])
    far_face_1 = last(cylinder_1[1])
    cylinder_1 = last(cylinder_1[2])

    # Remove duplicates but preserves tags for non-removed objects
    kernel.fragment(kernel.getEntities(), [])
    kernel.synchronize()

    boundaries = []
    for cylinder in [cylinder_0, cylinder_1]
        _, local_boundaries = gmsh.model.getAdjacencies(3, cylinder)
        local_idx = indexin(local_boundaries, boundaries)
        local_delete = []
        for (idx, boundary) in zip(local_idx, local_boundaries)
            if isnothing(idx)
                push!(boundaries, boundary)
            else
                push!(local_delete, idx)
            end
        end
        deleteat!(boundaries, sort(local_delete))
    end
    deleteat!(
        boundaries,
        sort(indexin([base_face_0, base_face_1, far_face_0, far_face_1], boundaries))
    )

    # Add physical groups
    cylinder_group =
        gmsh.model.addPhysicalGroup(3, [cylinder_0, cylinder_1], -1, "cylinder")
    boundary_group = gmsh.model.addPhysicalGroup(2, boundaries, -1, "boundaries")

    port1_group = gmsh.model.addPhysicalGroup(2, [base_face_0, base_face_1], -1, "port1")
    port2_group = gmsh.model.addPhysicalGroup(2, [far_face_0, far_face_1], -1, "port2")

    # Generate mesh
    gmsh.option.setNumber("Mesh.MinimumCurveNodes", 2)
    gmsh.option.setNumber("Mesh.MinimumCircleNodes", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    gmsh.model.mesh.generate(3)
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
        println("Port 1: ", port1_group)
        println("Port 2: ", port2_group)
        println()
    end

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
