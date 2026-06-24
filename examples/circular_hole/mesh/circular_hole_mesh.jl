# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Gmsh: gmsh

"""
    generate_circular_hole_mesh(;
        filename::AbstractString = "circular_hole.msh",
        L_outer::Real = 15.0,
        R::Real = 3.0,
        r::Real = 1.0,
        mesh_size_coarse::Real = 2.0,
        mesh_size_medium::Real = 0.4,
        mesh_size_fine::Real = 0.2,
        verbose::Integer = 5,
        gui::Bool = false
    )

Generate a mesh for a circular disk with a central hole example using Gmsh

# Arguments

  - filename - the filename to use for the generated mesh
  - L_outer - outer computational box size
  - R - radius of the circular disk
  - r - radius of the central circular hole
  - mesh_size_coarse - mesh size for box surfaces (coarse regions)
  - mesh_size_medium - mesh size for disk regions
  - mesh_size_fine - mesh size for hole regions (fine regions)
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_circular_hole_mesh(;
    filename::AbstractString="circular_hole.msh",
    L_outer::Real=15.0,
    R::Real=3.0,
    r::Real=1.0,
    mesh_size_coarse::Real=2.0,
    mesh_size_medium::Real=0.4,
    mesh_size_fine::Real=0.2,
    verbose::Integer=5,
    gui::Bool=false
)

    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    gmsh.model.add("circular_hole")

    # Set MSH file format version to 2.2
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

    # Center point
    cx, cy, cz = 0.5, 0.5, 0.5

    # Create cube vertices
    cube_pts = [
        gmsh.model.geo.addPoint(
            cx - L_outer / 2,
            cy - L_outer / 2,
            cz - L_outer / 2,
            mesh_size_coarse
        ),
        gmsh.model.geo.addPoint(
            cx + L_outer / 2,
            cy - L_outer / 2,
            cz - L_outer / 2,
            mesh_size_coarse
        ),
        gmsh.model.geo.addPoint(
            cx + L_outer / 2,
            cy + L_outer / 2,
            cz - L_outer / 2,
            mesh_size_coarse
        ),
        gmsh.model.geo.addPoint(
            cx - L_outer / 2,
            cy + L_outer / 2,
            cz - L_outer / 2,
            mesh_size_coarse
        ),
        gmsh.model.geo.addPoint(
            cx - L_outer / 2,
            cy - L_outer / 2,
            cz + L_outer / 2,
            mesh_size_coarse
        ),
        gmsh.model.geo.addPoint(
            cx + L_outer / 2,
            cy - L_outer / 2,
            cz + L_outer / 2,
            mesh_size_coarse
        ),
        gmsh.model.geo.addPoint(
            cx + L_outer / 2,
            cy + L_outer / 2,
            cz + L_outer / 2,
            mesh_size_coarse
        ),
        gmsh.model.geo.addPoint(
            cx - L_outer / 2,
            cy + L_outer / 2,
            cz + L_outer / 2,
            mesh_size_coarse
        )
    ]

    # Create disk points
    center_pt = gmsh.model.geo.addPoint(cx, cy, cz, mesh_size_medium)
    outer_pts = [
        gmsh.model.geo.addPoint(cx + R, cy, cz, mesh_size_medium),
        gmsh.model.geo.addPoint(cx, cy + R, cz, mesh_size_medium),
        gmsh.model.geo.addPoint(cx - R, cy, cz, mesh_size_medium),
        gmsh.model.geo.addPoint(cx, cy - R, cz, mesh_size_medium)
    ]
    inner_pts = [
        gmsh.model.geo.addPoint(cx + r, cy, cz, mesh_size_fine),
        gmsh.model.geo.addPoint(cx, cy + r, cz, mesh_size_fine),
        gmsh.model.geo.addPoint(cx - r, cy, cz, mesh_size_fine),
        gmsh.model.geo.addPoint(cx, cy - r, cz, mesh_size_fine)
    ]

    # Create cube edges
    cube_lines = [
        gmsh.model.geo.addLine(cube_pts[1], cube_pts[2]),
        gmsh.model.geo.addLine(cube_pts[2], cube_pts[3]),
        gmsh.model.geo.addLine(cube_pts[3], cube_pts[4]),
        gmsh.model.geo.addLine(cube_pts[4], cube_pts[1]),
        gmsh.model.geo.addLine(cube_pts[5], cube_pts[6]),
        gmsh.model.geo.addLine(cube_pts[6], cube_pts[7]),
        gmsh.model.geo.addLine(cube_pts[7], cube_pts[8]),
        gmsh.model.geo.addLine(cube_pts[8], cube_pts[5]),
        gmsh.model.geo.addLine(cube_pts[1], cube_pts[5]),
        gmsh.model.geo.addLine(cube_pts[2], cube_pts[6]),
        gmsh.model.geo.addLine(cube_pts[3], cube_pts[7]),
        gmsh.model.geo.addLine(cube_pts[4], cube_pts[8])
    ]

    # Create circles
    outer_circles = [
        gmsh.model.geo.addCircleArc(outer_pts[1], center_pt, outer_pts[2]),
        gmsh.model.geo.addCircleArc(outer_pts[2], center_pt, outer_pts[3]),
        gmsh.model.geo.addCircleArc(outer_pts[3], center_pt, outer_pts[4]),
        gmsh.model.geo.addCircleArc(outer_pts[4], center_pt, outer_pts[1])
    ]
    inner_circles = [
        gmsh.model.geo.addCircleArc(inner_pts[1], center_pt, inner_pts[2]),
        gmsh.model.geo.addCircleArc(inner_pts[2], center_pt, inner_pts[3]),
        gmsh.model.geo.addCircleArc(inner_pts[3], center_pt, inner_pts[4]),
        gmsh.model.geo.addCircleArc(inner_pts[4], center_pt, inner_pts[1])
    ]

    # Create curve loops and surfaces
    outer_loop = gmsh.model.geo.addCurveLoop(outer_circles)
    inner_loop = gmsh.model.geo.addCurveLoop(inner_circles)
    disk_surface = gmsh.model.geo.addPlaneSurface([outer_loop, inner_loop])
    hole_surface = gmsh.model.geo.addPlaneSurface([inner_loop])

    # Create cube surfaces
    cube_loops = [
        gmsh.model.geo.addCurveLoop([
            cube_lines[1],
            cube_lines[2],
            cube_lines[3],
            cube_lines[4]
        ]),
        gmsh.model.geo.addCurveLoop([
            cube_lines[5],
            cube_lines[6],
            cube_lines[7],
            cube_lines[8]
        ]),
        gmsh.model.geo.addCurveLoop([
            cube_lines[9],
            cube_lines[5],
            -cube_lines[10],
            -cube_lines[1]
        ]),
        gmsh.model.geo.addCurveLoop([
            cube_lines[10],
            cube_lines[6],
            -cube_lines[11],
            -cube_lines[2]
        ]),
        gmsh.model.geo.addCurveLoop([
            cube_lines[11],
            cube_lines[7],
            -cube_lines[12],
            -cube_lines[3]
        ]),
        gmsh.model.geo.addCurveLoop([
            cube_lines[12],
            cube_lines[8],
            -cube_lines[9],
            -cube_lines[4]
        ])
    ]
    cube_surfaces = [gmsh.model.geo.addPlaneSurface([loop]) for loop in cube_loops]

    # Create volume
    surface_loop =
        gmsh.model.geo.addSurfaceLoop([cube_surfaces..., disk_surface, hole_surface])
    volume = gmsh.model.geo.addVolume([surface_loop])

    # Synchronize geometry
    gmsh.model.geo.synchronize()

    # Add mesh size fields
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", inner_circles)

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", mesh_size_fine)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", mesh_size_medium)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.2)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 1.0)

    gmsh.model.mesh.field.add("Distance", 3)
    gmsh.model.mesh.field.setNumbers(3, "CurvesList", outer_circles)

    gmsh.model.mesh.field.add("Threshold", 4)
    gmsh.model.mesh.field.setNumber(4, "InField", 3)
    gmsh.model.mesh.field.setNumber(4, "SizeMin", mesh_size_medium)
    gmsh.model.mesh.field.setNumber(4, "SizeMax", mesh_size_coarse)
    gmsh.model.mesh.field.setNumber(4, "DistMin", 0.5)
    gmsh.model.mesh.field.setNumber(4, "DistMax", 2.5)

    # Add field to make cube surfaces coarse
    gmsh.model.mesh.field.add("Constant", 6)
    gmsh.model.mesh.field.setNumber(6, "VIn", mesh_size_coarse)
    gmsh.model.mesh.field.setNumber(6, "VOut", mesh_size_coarse)
    gmsh.model.mesh.field.setNumbers(6, "SurfacesList", cube_surfaces)

    gmsh.model.mesh.field.add("Min", 5)
    gmsh.model.mesh.field.setNumbers(5, "FieldsList", [2, 4])

    gmsh.model.mesh.field.add("Max", 7)
    gmsh.model.mesh.field.setNumbers(7, "FieldsList", [5, 6])
    gmsh.model.mesh.field.setAsBackgroundMesh(7)

    # Physical groups
    gmsh.model.addPhysicalGroup(3, [volume], 1, "domain")
    gmsh.model.addPhysicalGroup(2, [cube_surfaces[1]], 2, "cube_wall_m_xy")
    gmsh.model.addPhysicalGroup(2, [cube_surfaces[2]], 3, "cube_wall_p_xy")
    gmsh.model.addPhysicalGroup(2, [cube_surfaces[3]], 4, "cube_wall_m_xz")
    gmsh.model.addPhysicalGroup(2, [cube_surfaces[4]], 5, "cube_wall_p_yz")
    gmsh.model.addPhysicalGroup(2, [cube_surfaces[5]], 6, "cube_wall_p_xz")
    gmsh.model.addPhysicalGroup(2, [cube_surfaces[6]], 7, "cube_wall_m_yz")
    gmsh.model.addPhysicalGroup(2, [disk_surface], 8, "disk_surface")
    gmsh.model.addPhysicalGroup(2, [hole_surface], 9, "hole_surface")

    # Generate mesh
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.mesh.generate(3)
    gmsh.write(filename)

    # Print physical group information
    println("Generated mesh: ", filename)
    println("Domain: 1 (domain)")
    println("Box boundaries: 2-7")
    println("Disk surface: 8")
    println("Hole surface: 9")
    println()

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
