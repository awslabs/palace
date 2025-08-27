using Gmsh

# Initialize Gmsh
gmsh.initialize()

# Define the cube sizes
L_outer = 20.0
rect_width = 12.0
rect_height = 6.0
r = 1.0
d = 5.0

# Set mesh sizes for different regions
mesh_size_coarse = 2.0   # Cube surfaces (very coarse)
mesh_size_medium = 0.4   # Rectangle region
mesh_size_fine = 0.2     # Hole region

# Mesh algorithm and control
gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay algorithm for 3D mesh
gmsh.option.setNumber("Mesh.CharacteristicLengthFactor", 1.0)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.05)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 3.0)

# Define center point for the cube and the rectangle
cx = 0.5
cy = 0.5
cz = 0.5

# Define the points (vertices of the outer cube)
p1 = gmsh.model.geo.addPoint(cx-L_outer/2, cy-L_outer/2, cz-L_outer/2, mesh_size_coarse)
p2 = gmsh.model.geo.addPoint(cx+L_outer/2, cy-L_outer/2, cz-L_outer/2, mesh_size_coarse)
p3 = gmsh.model.geo.addPoint(cx+L_outer/2, cy+L_outer/2, cz-L_outer/2, mesh_size_coarse)
p4 = gmsh.model.geo.addPoint(cx-L_outer/2, cy+L_outer/2, cz-L_outer/2, mesh_size_coarse)
p5 = gmsh.model.geo.addPoint(cx-L_outer/2, cy-L_outer/2, cz+L_outer/2, mesh_size_coarse)
p6 = gmsh.model.geo.addPoint(cx+L_outer/2, cy-L_outer/2, cz+L_outer/2, mesh_size_coarse)
p7 = gmsh.model.geo.addPoint(cx+L_outer/2, cy+L_outer/2, cz+L_outer/2, mesh_size_coarse)
p8 = gmsh.model.geo.addPoint(cx-L_outer/2, cy+L_outer/2, cz+L_outer/2, mesh_size_coarse)

# Rectangle corners
p9 = gmsh.model.geo.addPoint(cx-rect_width/2, cy-rect_height/2, cz, mesh_size_medium)  # Bottom left
p10 = gmsh.model.geo.addPoint(cx+rect_width/2, cy-rect_height/2, cz, mesh_size_medium) # Bottom right
p11 = gmsh.model.geo.addPoint(cx+rect_width/2, cy+rect_height/2, cz, mesh_size_medium) # Top right
p12 = gmsh.model.geo.addPoint(cx-rect_width/2, cy+rect_height/2, cz, mesh_size_medium) # Top left

# Left hole center and points
p13 = gmsh.model.geo.addPoint(cx-d/2, cy, cz, mesh_size_fine) # Left hole center
p14 = gmsh.model.geo.addPoint(cx-d/2+r, cy, cz, mesh_size_fine) # Right point of left hole
p15 = gmsh.model.geo.addPoint(cx-d/2, cy+r, cz, mesh_size_fine) # Top point of left hole
p16 = gmsh.model.geo.addPoint(cx-d/2-r, cy, cz, mesh_size_fine) # Left point of left hole
p17 = gmsh.model.geo.addPoint(cx-d/2, cy-r, cz, mesh_size_fine) # Bottom point of left hole

# Right hole center and points
p18 = gmsh.model.geo.addPoint(cx+d/2, cy, cz, mesh_size_fine) # Right hole center
p19 = gmsh.model.geo.addPoint(cx+d/2+r, cy, cz, mesh_size_fine) # Right point of right hole
p20 = gmsh.model.geo.addPoint(cx+d/2, cy+r, cz, mesh_size_fine) # Top point of right hole
p21 = gmsh.model.geo.addPoint(cx+d/2-r, cy, cz, mesh_size_fine) # Left point of right hole
p22 = gmsh.model.geo.addPoint(cx+d/2, cy-r, cz, mesh_size_fine) # Bottom point of right hole

# Define the lines (edges of the outer cube)
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

# Rectangle edges
l13 = gmsh.model.geo.addLine(p9, p10)  # Bottom edge
l14 = gmsh.model.geo.addLine(p10, p11) # Right edge
l15 = gmsh.model.geo.addLine(p11, p12) # Top edge
l16 = gmsh.model.geo.addLine(p12, p9)  # Left edge

# Left hole circles
l17 = gmsh.model.geo.addCircleArc(p14, p13, p15) # Left hole quarter 1
l18 = gmsh.model.geo.addCircleArc(p15, p13, p16) # Left hole quarter 2
l19 = gmsh.model.geo.addCircleArc(p16, p13, p17) # Left hole quarter 3
l20 = gmsh.model.geo.addCircleArc(p17, p13, p14) # Left hole quarter 4

# Right hole circles
l21 = gmsh.model.geo.addCircleArc(p19, p18, p20) # Right hole quarter 1
l22 = gmsh.model.geo.addCircleArc(p20, p18, p21) # Right hole quarter 2
l23 = gmsh.model.geo.addCircleArc(p21, p18, p22) # Right hole quarter 3
l24 = gmsh.model.geo.addCircleArc(p22, p18, p19) # Right hole quarter 4

# Create curve loops
cl1 = gmsh.model.geo.addCurveLoop([l13, l14, l15, l16]) # Rectangle boundary
cl2 = gmsh.model.geo.addCurveLoop([l17, l18, l19, l20]) # Left hole
cl3 = gmsh.model.geo.addCurveLoop([l21, l22, l23, l24]) # Right hole

# Create rectangle surface with two holes
s1 = gmsh.model.geo.addPlaneSurface([cl1, cl2, cl3]) # Rectangle minus both holes
s2 = gmsh.model.geo.addPlaneSurface([cl2]) # Left hole surface
s3 = gmsh.model.geo.addPlaneSurface([cl3]) # Right hole surface

# Define the surfaces (faces of the outer cube)
cl4 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
s4 = gmsh.model.geo.addPlaneSurface([cl4])

cl5 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])
s5 = gmsh.model.geo.addPlaneSurface([cl5])

cl6 = gmsh.model.geo.addCurveLoop([l9, l5, -l10, -l1])
s6 = gmsh.model.geo.addPlaneSurface([cl6])

cl7 = gmsh.model.geo.addCurveLoop([l10, l6, -l11, -l2])
s7 = gmsh.model.geo.addPlaneSurface([cl7])

cl8 = gmsh.model.geo.addCurveLoop([l11, l7, -l12, -l3])
s8 = gmsh.model.geo.addPlaneSurface([cl8])

cl9 = gmsh.model.geo.addCurveLoop([l12, l8, -l9, -l4])
s9 = gmsh.model.geo.addPlaneSurface([cl9])

# Create cube volume with embedded rectangle surfaces
sl1 = gmsh.model.geo.addSurfaceLoop([s4, s5, s6, s7, s8, s9, s1, s2, s3]) # Include rectangle surfaces
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
gmsh.model.mesh.field.setNumbers(field3, "CurvesList", [l21, l22, l23, l24])

field4 = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(field4, "InField", field3)
gmsh.model.mesh.field.setNumber(field4, "SizeMin", mesh_size_fine)
gmsh.model.mesh.field.setNumber(field4, "SizeMax", mesh_size_medium)
gmsh.model.mesh.field.setNumber(field4, "DistMin", 0.2)
gmsh.model.mesh.field.setNumber(field4, "DistMax", 1.0)

# Field for rectangle edges
field5 = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(field5, "CurvesList", [l13, l14, l15, l16])

field6 = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(field6, "InField", field5)
gmsh.model.mesh.field.setNumber(field6, "SizeMin", mesh_size_medium)
gmsh.model.mesh.field.setNumber(field6, "SizeMax", mesh_size_coarse)
gmsh.model.mesh.field.setNumber(field6, "DistMin", 0.5)
gmsh.model.mesh.field.setNumber(field6, "DistMax", 2.5)

# Field to make cube surfaces coarse
field7 = gmsh.model.mesh.field.add("Constant")
gmsh.model.mesh.field.setNumber(field7, "VIn", mesh_size_coarse)
gmsh.model.mesh.field.setNumber(field7, "VOut", mesh_size_coarse)
gmsh.model.mesh.field.setNumbers(field7, "SurfacesList", [s4, s5, s6, s7, s8, s9])

# Combine fields
field8 = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(field8, "FieldsList", [field2, field4, field6])

field9 = gmsh.model.mesh.field.add("Max")
gmsh.model.mesh.field.setNumbers(field9, "FieldsList", [field8, field7])

gmsh.model.mesh.field.setAsBackgroundMesh(field9)

# Physical groups
gmsh.model.addPhysicalGroup(3, [v1], 1)
gmsh.model.setPhysicalName(3, 1, "domain")

gmsh.model.addPhysicalGroup(2, [s4], 2)
gmsh.model.setPhysicalName(2, 2, "cube_wall_m_xy")

gmsh.model.addPhysicalGroup(2, [s5], 3)
gmsh.model.setPhysicalName(2, 3, "cube_wall_p_xy")

gmsh.model.addPhysicalGroup(2, [s6], 4)
gmsh.model.setPhysicalName(2, 4, "cube_wall_m_xz")

gmsh.model.addPhysicalGroup(2, [s7], 5)
gmsh.model.setPhysicalName(2, 5, "cube_wall_p_yz")

gmsh.model.addPhysicalGroup(2, [s8], 6)
gmsh.model.setPhysicalName(2, 6, "cube_wall_p_xz")

gmsh.model.addPhysicalGroup(2, [s9], 7)
gmsh.model.setPhysicalName(2, 7, "cube_Wall_m_yz")

gmsh.model.addPhysicalGroup(2, [s1], 8)
gmsh.model.setPhysicalName(2, 8, "rectangle_surface")

gmsh.model.addPhysicalGroup(2, [s2], 9)
gmsh.model.setPhysicalName(2, 9, "left_hole_surface")

gmsh.model.addPhysicalGroup(2, [s3], 10)
gmsh.model.setPhysicalName(2, 10, "right_hole_surface")

# Generate mesh
gmsh.model.mesh.generate(3)

# Write mesh file
gmsh.write("sheet_w_two_holes.msh")

# Finalize Gmsh
gmsh.finalize()
