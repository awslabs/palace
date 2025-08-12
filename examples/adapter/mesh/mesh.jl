# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Gmsh: gmsh

"""
    generate_adapter_mesh(;
        filename::AbstractString,
        refinement::Integer  = 0,
        order::Integer       = 2,
        ellipse_a::Real      = 1.5,
        ellipse_b::Real      = 0.9,
        rect_a::Real         = 2.3,
        rect_b::Real         = 1.0,
        adapter_length::Real = 14.0,
        verbose::Integer     = 5,
        gui::Bool            = false
    )

Generate a mesh for the adapter waveguide example using Gmsh

# Arguments

  - filename - the filename to use for the generated mesh
  - refinement - measure of how many elements to include, 0 is least
  - order - the polynomial order of the approximation, minimum 1
  - ellipse_a - x radius of ellipse
  - ellipse_b - y radius of ellipse
  - rect_a - x dimension of rectangle
  - rect_b - y dimension of rectangle
  - adapter_length - the length of the adapter waveguide
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_adapter_mesh(;
    filename::AbstractString,
    refinement::Integer  = 0,
    order::Integer       = 2,
    ellipse_a::Real      = 1.5,
    ellipse_b::Real      = 0.9,
    rect_a::Real         = 2.3,
    rect_b::Real         = 1.0,
    adapter_length::Real = 14.0,
    verbose::Integer     = 5,
    gui::Bool            = false
)
    @assert refinement >= 0
    @assert order > 0

    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "adapter" in gmsh.model.list()
        gmsh.model.setCurrent("adapter")
        gmsh.model.remove()
    end
    gmsh.model.add("adapter")

    # Geometry
    # Create ellipse points and arcs
    # Split the ellipse into 8 segments
    xp = rect_a / 2.0
    yp = ellipse_b * sqrt(1.0 - xp^2 / ellipse_a^2)

    ellipse_points = [
        kernel.addPoint(ellipse_a, 0, 0),
        kernel.addPoint(xp, yp, 0),
        kernel.addPoint(0, ellipse_b, 0),
        kernel.addPoint(-xp, yp, 0),
        kernel.addPoint(-ellipse_a, 0, 0),
        kernel.addPoint(-xp, -yp, 0),
        kernel.addPoint(0, -ellipse_b, 0),
        kernel.addPoint(xp, -yp, 0)
    ]
    center = kernel.addPoint(0, 0, 0)

    ellipse_arcs = [
        kernel.addEllipseArc(
            ellipse_points[1],
            center,
            ellipse_points[1],
            ellipse_points[2]
        ),
        kernel.addEllipseArc(
            ellipse_points[2],
            center,
            ellipse_points[3],
            ellipse_points[3]
        ),
        kernel.addEllipseArc(
            ellipse_points[3],
            center,
            ellipse_points[3],
            ellipse_points[4]
        ),
        kernel.addEllipseArc(
            ellipse_points[4],
            center,
            ellipse_points[5],
            ellipse_points[5]
        ),
        kernel.addEllipseArc(
            ellipse_points[5],
            center,
            ellipse_points[5],
            ellipse_points[6]
        ),
        kernel.addEllipseArc(
            ellipse_points[6],
            center,
            ellipse_points[7],
            ellipse_points[7]
        ),
        kernel.addEllipseArc(
            ellipse_points[7],
            center,
            ellipse_points[7],
            ellipse_points[8]
        ),
        kernel.addEllipseArc(
            ellipse_points[8],
            center,
            ellipse_points[1],
            ellipse_points[1]
        )
    ]

    # Create the ellipse surface
    ellipse_loop = kernel.addCurveLoop(ellipse_arcs)
    ellipse = kernel.addPlaneSurface([ellipse_loop])

    # Create rectangle points and lines
    # Split rectangle into 8 segments
    rect_points = [
        kernel.addPoint(0.5 * rect_a, 0, adapter_length),
        kernel.addPoint(0.5 * rect_a, 0.5 * rect_b, adapter_length),
        kernel.addPoint(0, 0.5 * rect_b, adapter_length),
        kernel.addPoint(-0.5 * rect_a, 0.5 * rect_b, adapter_length),
        kernel.addPoint(-0.5 * rect_a, 0, adapter_length),
        kernel.addPoint(-0.5 * rect_a, -0.5 * rect_b, adapter_length),
        kernel.addPoint(0, -0.5 * rect_b, adapter_length),
        kernel.addPoint(0.5 * rect_a, -0.5 * rect_b, adapter_length)
    ]

    rect_lines = [
        kernel.addLine(rect_points[1], rect_points[2]),
        kernel.addLine(rect_points[2], rect_points[3]),
        kernel.addLine(rect_points[3], rect_points[4]),
        kernel.addLine(rect_points[4], rect_points[5]),
        kernel.addLine(rect_points[5], rect_points[6]),
        kernel.addLine(rect_points[6], rect_points[7]),
        kernel.addLine(rect_points[7], rect_points[8]),
        kernel.addLine(rect_points[8], rect_points[1])
    ]

    # Create rectangle surface
    rect_loop = kernel.addCurveLoop(rect_lines)
    rect = kernel.addPlaneSurface([rect_loop])

    # Create volume
    volume = kernel.addThruSections([ellipse, rect], -1, true, false, -1, "C0")
    volume = last.(volume)
    volume_idx = first(volume)

    # Fragment volume
    kernel.fragment([(3, volume_idx)], [(2, rect), (2, ellipse)])

    kernel.synchronize()

    # Add physical groups
    adapter_group = gmsh.model.addPhysicalGroup(3, volume, -1, "adapter")

    bottom = last.(
        gmsh.model.getEntitiesInBoundingBox(
            -2 * ellipse_a,
            -2 * ellipse_b,
            -0.1 * adapter_length,
            2 * ellipse_a,
            2 * ellipse_b,
            0.5 * adapter_length,
            2
        )
    )
    top = last.(
        gmsh.model.getEntitiesInBoundingBox(
            -0.6 * rect_a,
            -0.6 * rect_b,
            0.5 * adapter_length,
            0.6 * rect_a,
            0.6 * rect_b,
            1.1 * adapter_length,
            2
        )
    )

    exterior = []
    for domain in volume
        _, domain_boundaries = gmsh.model.getAdjacencies(3, domain)
        for boundary in domain_boundaries
            if boundary in top || boundary in bottom
                continue
            end
            idx = first(indexin([boundary], exterior))
            if isnothing(idx)
                push!(exterior, boundary)
            else
                deleteat!(exterior, idx)
            end
        end
    end
    top_group = gmsh.model.addPhysicalGroup(2, top, -1, "top")
    bottom_group = gmsh.model.addPhysicalGroup(2, bottom, -1, "bottom")
    exterior_group = gmsh.model.addPhysicalGroup(2, exterior, -1, "exterior")

    # Generate mesh
    gmsh.option.setNumber(
        "Mesh.MeshSizeMax",
        minimum([ellipse_a, ellipse_b, rect_a/2, rect_b/2, adapter_length])
    )
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)

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
        println("Adapter: ", adapter_group)
        println("Exterior boundaries: ", exterior_group)
        println("Top boundaries: ", top_group)
        println("Bottom boundaries: ", bottom_group)
    end

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
