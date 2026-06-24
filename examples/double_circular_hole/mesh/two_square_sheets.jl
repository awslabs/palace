# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Gmsh: gmsh

"""
    generate_two_square_sheets_mesh(;
        filename::AbstractString = "two_square_sheets2.msh",
        L_outer::Real = 20.0,
        sheet_size::Real = 4.0,
        hole_radius::Real = 0.8,
        sheet_separation::Real = 6.0,
        mesh_size_coarse::Real = 3.0,
        mesh_size_medium::Real = 0.5,
        mesh_size_fine::Real = 0.2,
        verbose::Integer = 5,
        gui::Bool = false
    )

Generate a mesh for the two square sheets with circular holes example using Gmsh

# Arguments

  - filename - the filename to use for the generated mesh
  - L_outer - outer computational box size
  - sheet_size - side length of each square sheet
  - hole_radius - radius of circular holes in the sheets
  - sheet_separation - distance between sheet centers
  - mesh_size_coarse - mesh size for box surfaces (coarse regions)
  - mesh_size_medium - mesh size for sheet regions
  - mesh_size_fine - mesh size for hole regions (fine regions)
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_two_square_sheets_mesh(;
    filename::AbstractString="two_square_sheets.msh",
    L_outer::Real=20.0,
    sheet_size::Real=4.0,
    hole_radius::Real=0.8,
    sheet_separation::Real=6.0,
    mesh_size_coarse::Real=3.0,
    mesh_size_medium::Real=0.5,
    mesh_size_fine::Real=0.2,
    verbose::Integer=5,
    gui::Bool=false
)

    # Initialize Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Mesh algorithm and control
    gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", 1.0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.05)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 3.0)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    # Define center point for the box and the sheets
    cx = 0.5
    cy = 0.5
    cz = 0.5

    # Define the points (vertices of the outer box)
    p1 = gmsh.model.geo.addPoint(
        cx - L_outer / 2,
        cy - L_outer / 2,
        cz - L_outer / 2,
        mesh_size_coarse
    )
    p2 = gmsh.model.geo.addPoint(
        cx + L_outer / 2,
        cy - L_outer / 2,
        cz - L_outer / 2,
        mesh_size_coarse
    )
    p3 = gmsh.model.geo.addPoint(
        cx + L_outer / 2,
        cy + L_outer / 2,
        cz - L_outer / 2,
        mesh_size_coarse
    )
    p4 = gmsh.model.geo.addPoint(
        cx - L_outer / 2,
        cy + L_outer / 2,
        cz - L_outer / 2,
        mesh_size_coarse
    )
    p5 = gmsh.model.geo.addPoint(
        cx - L_outer / 2,
        cy - L_outer / 2,
        cz + L_outer / 2,
        mesh_size_coarse
    )
    p6 = gmsh.model.geo.addPoint(
        cx + L_outer / 2,
        cy - L_outer / 2,
        cz + L_outer / 2,
        mesh_size_coarse
    )
    p7 = gmsh.model.geo.addPoint(
        cx + L_outer / 2,
        cy + L_outer / 2,
        cz + L_outer / 2,
        mesh_size_coarse
    )
    p8 = gmsh.model.geo.addPoint(
        cx - L_outer / 2,
        cy + L_outer / 2,
        cz + L_outer / 2,
        mesh_size_coarse
    )

    # Left sheet corner points
    p9 = gmsh.model.geo.addPoint(
        cx - sheet_separation / 2 - sheet_size / 2,
        cy - sheet_size / 2,
        cz,
        mesh_size_medium
    )  # Bottom left
    p10 = gmsh.model.geo.addPoint(
        cx - sheet_separation / 2 + sheet_size / 2,
        cy - sheet_size / 2,
        cz,
        mesh_size_medium
    ) # Bottom right
    p11 = gmsh.model.geo.addPoint(
        cx - sheet_separation / 2 + sheet_size / 2,
        cy + sheet_size / 2,
        cz,
        mesh_size_medium
    ) # Top right
    p12 = gmsh.model.geo.addPoint(
        cx - sheet_separation / 2 - sheet_size / 2,
        cy + sheet_size / 2,
        cz,
        mesh_size_medium
    ) # Top left

    # Left hole center and points
    p13 = gmsh.model.geo.addPoint(cx - sheet_separation / 2, cy, cz, mesh_size_fine) # Left hole center
    p14 = gmsh.model.geo.addPoint(
        cx - sheet_separation / 2 + hole_radius,
        cy,
        cz,
        mesh_size_fine
    ) # Right point of left hole
    p15 = gmsh.model.geo.addPoint(
        cx - sheet_separation / 2,
        cy + hole_radius,
        cz,
        mesh_size_fine
    ) # Top point of left hole
    p16 = gmsh.model.geo.addPoint(
        cx - sheet_separation / 2 - hole_radius,
        cy,
        cz,
        mesh_size_fine
    ) # Left point of left hole
    p17 = gmsh.model.geo.addPoint(
        cx - sheet_separation / 2,
        cy - hole_radius,
        cz,
        mesh_size_fine
    ) # Bottom point of left hole

    # Right sheet corner points
    p18 = gmsh.model.geo.addPoint(
        cx + sheet_separation / 2 - sheet_size / 2,
        cy - sheet_size / 2,
        cz,
        mesh_size_medium
    ) # Bottom left
    p19 = gmsh.model.geo.addPoint(
        cx + sheet_separation / 2 + sheet_size / 2,
        cy - sheet_size / 2,
        cz,
        mesh_size_medium
    ) # Bottom right
    p20 = gmsh.model.geo.addPoint(
        cx + sheet_separation / 2 + sheet_size / 2,
        cy + sheet_size / 2,
        cz,
        mesh_size_medium
    ) # Top right
    p21 = gmsh.model.geo.addPoint(
        cx + sheet_separation / 2 - sheet_size / 2,
        cy + sheet_size / 2,
        cz,
        mesh_size_medium
    ) # Top left

    # Right hole center and points
    p22 = gmsh.model.geo.addPoint(cx + sheet_separation / 2, cy, cz, mesh_size_fine) # Right hole center
    p23 = gmsh.model.geo.addPoint(
        cx + sheet_separation / 2 + hole_radius,
        cy,
        cz,
        mesh_size_fine
    ) # Right point of right hole
    p24 = gmsh.model.geo.addPoint(
        cx + sheet_separation / 2,
        cy + hole_radius,
        cz,
        mesh_size_fine
    ) # Top point of right hole
    p25 = gmsh.model.geo.addPoint(
        cx + sheet_separation / 2 - hole_radius,
        cy,
        cz,
        mesh_size_fine
    ) # Left point of right hole
    p26 = gmsh.model.geo.addPoint(
        cx + sheet_separation / 2,
        cy - hole_radius,
        cz,
        mesh_size_fine
    ) # Bottom point of right hole

    # Define the lines (edges of the outer box)
    l1 = gmsh.model.geo.addLine(p1, p2)
    l2 = gmsh.model.geo.addLine(p2, p3)
    l3 = gmsh.model.geo.addLine(p3, p4)
    l4 = gmsh.model.geo.addLine(p4, p1)
    l5 = gmsh.model.geo.addLine(p5, p6)
    l6 = gmsh.model.geo.addLine(p6, p7)
    l7 = gmsh.model.geo.addLine(p7, p8)
    l8 = gmsh.model.geo.addLine(p8, p5)
    l9 = gmsh.model.geo.addLine(p1, p5)
    l10 = gmsh.model.geo.addLine(p2, p6)
    l11 = gmsh.model.geo.addLine(p3, p7)
    l12 = gmsh.model.geo.addLine(p4, p8)

    # Left sheet edges
    l13 = gmsh.model.geo.addLine(p9, p10)  # Bottom edge
    l14 = gmsh.model.geo.addLine(p10, p11) # Right edge  
    l15 = gmsh.model.geo.addLine(p11, p12) # Top edge
    l16 = gmsh.model.geo.addLine(p12, p9)  # Left edge

    # Left hole circles
    l17 = gmsh.model.geo.addCircleArc(p14, p13, p15) # Left hole quarter 1
    l18 = gmsh.model.geo.addCircleArc(p15, p13, p16) # Left hole quarter 2
    l19 = gmsh.model.geo.addCircleArc(p16, p13, p17) # Left hole quarter 3
    l20 = gmsh.model.geo.addCircleArc(p17, p13, p14) # Left hole quarter 4

    # Right sheet edges
    l21 = gmsh.model.geo.addLine(p18, p19) # Bottom edge
    l22 = gmsh.model.geo.addLine(p19, p20) # Right edge
    l23 = gmsh.model.geo.addLine(p20, p21) # Top edge
    l24 = gmsh.model.geo.addLine(p21, p18) # Left edge

    # Right hole circles
    l25 = gmsh.model.geo.addCircleArc(p23, p22, p24) # Right hole quarter 1
    l26 = gmsh.model.geo.addCircleArc(p24, p22, p25) # Right hole quarter 2
    l27 = gmsh.model.geo.addCircleArc(p25, p22, p26) # Right hole quarter 3
    l28 = gmsh.model.geo.addCircleArc(p26, p22, p23) # Right hole quarter 4

    # Create curve loops
    cl1 = gmsh.model.geo.addCurveLoop([l13, l14, l15, l16]) # Left sheet boundary
    cl2 = gmsh.model.geo.addCurveLoop([l17, l18, l19, l20]) # Left hole
    cl3 = gmsh.model.geo.addCurveLoop([l21, l22, l23, l24]) # Right sheet boundary
    cl4 = gmsh.model.geo.addCurveLoop([l25, l26, l27, l28]) # Right hole

    # Create sheet surfaces with holes
    s1 = gmsh.model.geo.addPlaneSurface([cl1, cl2]) # Left sheet minus left hole
    s2 = gmsh.model.geo.addPlaneSurface([cl2]) # Left hole surface
    s3 = gmsh.model.geo.addPlaneSurface([cl3, cl4]) # Right sheet minus right hole
    s4 = gmsh.model.geo.addPlaneSurface([cl4]) # Right hole surface

    # Define the surfaces (faces of the outer box)
    cl5 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
    s5 = gmsh.model.geo.addPlaneSurface([cl5])

    cl6 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])
    s6 = gmsh.model.geo.addPlaneSurface([cl6])

    cl7 = gmsh.model.geo.addCurveLoop([l9, l5, -l10, -l1])
    s7 = gmsh.model.geo.addPlaneSurface([cl7])

    cl8 = gmsh.model.geo.addCurveLoop([l10, l6, -l11, -l2])
    s8 = gmsh.model.geo.addPlaneSurface([cl8])

    cl9 = gmsh.model.geo.addCurveLoop([l11, l7, -l12, -l3])
    s9 = gmsh.model.geo.addPlaneSurface([cl9])

    cl10 = gmsh.model.geo.addCurveLoop([l12, l8, -l9, -l4])
    s10 = gmsh.model.geo.addPlaneSurface([cl10])

    # Create box volume with embedded sheet surfaces
    sl1 = gmsh.model.geo.addSurfaceLoop([s5, s6, s7, s8, s9, s10, s1, s2, s3, s4]) # Include sheet surfaces
    v1 = gmsh.model.geo.addVolume([sl1])

    # Synchronize the model
    gmsh.model.geo.synchronize()

    # Mesh size fields for better control
    # Field for left hole edges
    field1 = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(field1, "CurvesList", [l17, l18, l19, l20])

    field2 = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(field2, "InField", field1)
    gmsh.model.mesh.field.setNumber(field2, "SizeMin", mesh_size_fine)
    gmsh.model.mesh.field.setNumber(field2, "SizeMax", mesh_size_medium)
    gmsh.model.mesh.field.setNumber(field2, "DistMin", 0.2)
    gmsh.model.mesh.field.setNumber(field2, "DistMax", 1.0)

    # Field for right hole edges
    field3 = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(field3, "CurvesList", [l25, l26, l27, l28])

    field4 = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(field4, "InField", field3)
    gmsh.model.mesh.field.setNumber(field4, "SizeMin", mesh_size_fine)
    gmsh.model.mesh.field.setNumber(field4, "SizeMax", mesh_size_medium)
    gmsh.model.mesh.field.setNumber(field4, "DistMin", 0.2)
    gmsh.model.mesh.field.setNumber(field4, "DistMax", 1.0)

    # Field for sheet edges
    field5 = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(
        field5,
        "CurvesList",
        [l13, l14, l15, l16, l21, l22, l23, l24]
    )

    field6 = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(field6, "InField", field5)
    gmsh.model.mesh.field.setNumber(field6, "SizeMin", mesh_size_medium)
    gmsh.model.mesh.field.setNumber(field6, "SizeMax", mesh_size_coarse)
    gmsh.model.mesh.field.setNumber(field6, "DistMin", 0.5)
    gmsh.model.mesh.field.setNumber(field6, "DistMax", 2.5)

    # Field to make box surfaces coarse
    field7 = gmsh.model.mesh.field.add("Constant")
    gmsh.model.mesh.field.setNumber(field7, "VIn", mesh_size_coarse)
    gmsh.model.mesh.field.setNumber(field7, "VOut", mesh_size_coarse)
    gmsh.model.mesh.field.setNumbers(field7, "SurfacesList", [s5, s6, s7, s8, s9, s10])

    # Combine fields
    field8 = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(field8, "FieldsList", [field2, field4, field6])

    field9 = gmsh.model.mesh.field.add("Max")
    gmsh.model.mesh.field.setNumbers(field9, "FieldsList", [field8, field7])

    gmsh.model.mesh.field.setAsBackgroundMesh(field9)

    # Physical groups
    gmsh.model.addPhysicalGroup(3, [v1], 1)
    gmsh.model.setPhysicalName(3, 1, "domain")

    gmsh.model.addPhysicalGroup(2, [s5], 2)
    gmsh.model.setPhysicalName(2, 2, "box_wall_m_xy")

    gmsh.model.addPhysicalGroup(2, [s6], 3)
    gmsh.model.setPhysicalName(2, 3, "box_wall_p_xy")

    gmsh.model.addPhysicalGroup(2, [s7], 4)
    gmsh.model.setPhysicalName(2, 4, "box_wall_m_xz")

    gmsh.model.addPhysicalGroup(2, [s8], 5)
    gmsh.model.setPhysicalName(2, 5, "box_wall_p_yz")

    gmsh.model.addPhysicalGroup(2, [s9], 6)
    gmsh.model.setPhysicalName(2, 6, "box_wall_p_xz")

    gmsh.model.addPhysicalGroup(2, [s10], 7)
    gmsh.model.setPhysicalName(2, 7, "box_wall_m_yz")

    gmsh.model.addPhysicalGroup(2, [s1], 8)
    gmsh.model.setPhysicalName(2, 8, "left_sheet_surface")

    gmsh.model.addPhysicalGroup(2, [s2], 9)
    gmsh.model.setPhysicalName(2, 9, "left_hole_surface")

    gmsh.model.addPhysicalGroup(2, [s3], 10)
    gmsh.model.setPhysicalName(2, 10, "right_sheet_surface")

    gmsh.model.addPhysicalGroup(2, [s4], 11)
    gmsh.model.setPhysicalName(2, 11, "right_hole_surface")
    # Generate mesh
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(2)

    # Write mesh file
    gmsh.write(filename)

    # Print physical group information
    println("Generated mesh: ", filename)
    println("Domain: 1 (domain)")
    println("Box boundaries: 2-7")
    println("Left sheet surface: 8")
    println("Left hole surface: 9")
    println("Right sheet surface: 10")
    println("Right hole surface: 11")
    println()

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
