# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Gmsh: gmsh

"""
    generate_rect_hole_mesh(;
        filename = "rect_hole.msh",
        dx = 2.0, dy = 5.0,      # hole side lengths
        W = 3.0,                 # film border width (outer film = (dx+2W) x (dy+2W))
        L_outer = 30.0,          # computational box size
        mesh_size_coarse = 3.0,
        mesh_size_medium = 0.5,
        mesh_size_fine = 0.2,
        verbose = 5, gui = false)

Rectangular superconducting film with a central rectangular hole (Khapaev 1997,
Supercond. Sci. Technol. 10 389, table 1 / figure 3). Attribute scheme matches
circular_hole: domain=1, box walls=2-7, film=8, hole=9.
"""
function generate_rect_hole_mesh(;
    filename::AbstractString="rect_hole.msh",
    dx::Real=2.0,
    dy::Real=5.0,
    W::Real=3.0,
    L_outer::Real=30.0,
    mesh_size_coarse::Real=3.0,
    mesh_size_medium::Real=0.5,
    mesh_size_fine::Real=0.2,
    verbose::Integer=5,
    gui::Bool=false
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    gmsh.model.add("rect_hole")
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

    cx, cy, cz = 0.5, 0.5, 0.5
    hx, hy = dx / 2, dy / 2            # hole half-sides
    fx, fy = dx / 2 + W, dy / 2 + W    # film half-sides

    cube_pts = [
        gmsh.model.geo.addPoint(cx - L_outer / 2, cy - L_outer / 2, cz - L_outer / 2, mesh_size_coarse),
        gmsh.model.geo.addPoint(cx + L_outer / 2, cy - L_outer / 2, cz - L_outer / 2, mesh_size_coarse),
        gmsh.model.geo.addPoint(cx + L_outer / 2, cy + L_outer / 2, cz - L_outer / 2, mesh_size_coarse),
        gmsh.model.geo.addPoint(cx - L_outer / 2, cy + L_outer / 2, cz - L_outer / 2, mesh_size_coarse),
        gmsh.model.geo.addPoint(cx - L_outer / 2, cy - L_outer / 2, cz + L_outer / 2, mesh_size_coarse),
        gmsh.model.geo.addPoint(cx + L_outer / 2, cy - L_outer / 2, cz + L_outer / 2, mesh_size_coarse),
        gmsh.model.geo.addPoint(cx + L_outer / 2, cy + L_outer / 2, cz + L_outer / 2, mesh_size_coarse),
        gmsh.model.geo.addPoint(cx - L_outer / 2, cy + L_outer / 2, cz + L_outer / 2, mesh_size_coarse)
    ]

    outer_pts = [
        gmsh.model.geo.addPoint(cx - fx, cy - fy, cz, mesh_size_medium),
        gmsh.model.geo.addPoint(cx + fx, cy - fy, cz, mesh_size_medium),
        gmsh.model.geo.addPoint(cx + fx, cy + fy, cz, mesh_size_medium),
        gmsh.model.geo.addPoint(cx - fx, cy + fy, cz, mesh_size_medium)
    ]
    inner_pts = [
        gmsh.model.geo.addPoint(cx - hx, cy - hy, cz, mesh_size_fine),
        gmsh.model.geo.addPoint(cx + hx, cy - hy, cz, mesh_size_fine),
        gmsh.model.geo.addPoint(cx + hx, cy + hy, cz, mesh_size_fine),
        gmsh.model.geo.addPoint(cx - hx, cy + hy, cz, mesh_size_fine)
    ]

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
    outer_lines = [
        gmsh.model.geo.addLine(outer_pts[1], outer_pts[2]),
        gmsh.model.geo.addLine(outer_pts[2], outer_pts[3]),
        gmsh.model.geo.addLine(outer_pts[3], outer_pts[4]),
        gmsh.model.geo.addLine(outer_pts[4], outer_pts[1])
    ]
    inner_lines = [
        gmsh.model.geo.addLine(inner_pts[1], inner_pts[2]),
        gmsh.model.geo.addLine(inner_pts[2], inner_pts[3]),
        gmsh.model.geo.addLine(inner_pts[3], inner_pts[4]),
        gmsh.model.geo.addLine(inner_pts[4], inner_pts[1])
    ]

    outer_loop = gmsh.model.geo.addCurveLoop(outer_lines)
    inner_loop = gmsh.model.geo.addCurveLoop(inner_lines)
    disk_surface = gmsh.model.geo.addPlaneSurface([outer_loop, inner_loop])
    hole_surface = gmsh.model.geo.addPlaneSurface([inner_loop])

    cube_loops = [
        gmsh.model.geo.addCurveLoop([cube_lines[1], cube_lines[2], cube_lines[3], cube_lines[4]]),
        gmsh.model.geo.addCurveLoop([cube_lines[5], cube_lines[6], cube_lines[7], cube_lines[8]]),
        gmsh.model.geo.addCurveLoop([cube_lines[9], cube_lines[5], -cube_lines[10], -cube_lines[1]]),
        gmsh.model.geo.addCurveLoop([cube_lines[10], cube_lines[6], -cube_lines[11], -cube_lines[2]]),
        gmsh.model.geo.addCurveLoop([cube_lines[11], cube_lines[7], -cube_lines[12], -cube_lines[3]]),
        gmsh.model.geo.addCurveLoop([cube_lines[12], cube_lines[8], -cube_lines[9], -cube_lines[4]])
    ]
    cube_surfaces = [gmsh.model.geo.addPlaneSurface([loop]) for loop in cube_loops]

    surface_loop =
        gmsh.model.geo.addSurfaceLoop([cube_surfaces..., disk_surface, hole_surface])
    volume = gmsh.model.geo.addVolume([surface_loop])
    gmsh.model.geo.synchronize()

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", inner_lines)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", mesh_size_fine)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", mesh_size_medium)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 1.5)

    gmsh.model.mesh.field.add("Distance", 3)
    gmsh.model.mesh.field.setNumbers(3, "CurvesList", outer_lines)
    gmsh.model.mesh.field.add("Threshold", 4)
    gmsh.model.mesh.field.setNumber(4, "InField", 3)
    gmsh.model.mesh.field.setNumber(4, "SizeMin", mesh_size_medium)
    gmsh.model.mesh.field.setNumber(4, "SizeMax", mesh_size_coarse)
    gmsh.model.mesh.field.setNumber(4, "DistMin", 0.5)
    gmsh.model.mesh.field.setNumber(4, "DistMax", 3.0)

    gmsh.model.mesh.field.add("Constant", 6)
    gmsh.model.mesh.field.setNumber(6, "VIn", mesh_size_coarse)
    gmsh.model.mesh.field.setNumber(6, "VOut", mesh_size_coarse)
    gmsh.model.mesh.field.setNumbers(6, "SurfacesList", cube_surfaces)

    gmsh.model.mesh.field.add("Min", 5)
    gmsh.model.mesh.field.setNumbers(5, "FieldsList", [2, 4])
    gmsh.model.mesh.field.add("Max", 7)
    gmsh.model.mesh.field.setNumbers(7, "FieldsList", [5, 6])
    gmsh.model.mesh.field.setAsBackgroundMesh(7)

    gmsh.model.addPhysicalGroup(3, [volume], 1, "domain")
    gmsh.model.addPhysicalGroup(2, [cube_surfaces[1]], 2, "cube_wall_m_xy")
    gmsh.model.addPhysicalGroup(2, [cube_surfaces[2]], 3, "cube_wall_p_xy")
    gmsh.model.addPhysicalGroup(2, [cube_surfaces[3]], 4, "cube_wall_m_xz")
    gmsh.model.addPhysicalGroup(2, [cube_surfaces[4]], 5, "cube_wall_p_yz")
    gmsh.model.addPhysicalGroup(2, [cube_surfaces[5]], 6, "cube_wall_p_xz")
    gmsh.model.addPhysicalGroup(2, [cube_surfaces[6]], 7, "cube_wall_m_yz")
    gmsh.model.addPhysicalGroup(2, [disk_surface], 8, "film_surface")
    gmsh.model.addPhysicalGroup(2, [hole_surface], 9, "hole_surface")

    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.model.mesh.generate(3)
    gmsh.write(filename)

    println("Generated mesh: ", filename)
    println("  hole = ", dx, " x ", dy, ", border W = ", W)
    println("  outer film = ", dx + 2W, " x ", dy + 2W, ", box = ", L_outer)
    println("Domain: 1  Box walls: 2-7  Film: 8  Hole: 9")

    gui && gmsh.fltk.run()
    return gmsh.finalize()
end
