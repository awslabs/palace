// Define the cube sizes
L_outer = 20.0;
rect_width = 12.0;
rect_height = 6.0;
r = 1.0;
d = 5.0;

// Set mesh sizes for different regions
mesh_size_coarse = 2.0;  // Cube surfaces (very coarse)
mesh_size_medium = 0.4;  // Rectangle region
mesh_size_fine = 0.2;    // Hole region

// Mesh algorithm and control
Mesh.Algorithm3D = 1; // Delaunay algorithm for 3D mesh
Mesh.CharacteristicLengthFactor = 1.0;
Mesh.MshFileVersion = 2.2;
Mesh.CharacteristicLengthMin = 0.05;
Mesh.CharacteristicLengthMax = 3.0;

// Define center point for the cube and the rectangle
cx = 0.5;
cy = 0.5;
cz = 0.5;

// Define the points (vertices of the outer cube)
Point(1) = {cx-L_outer/2, cy-L_outer/2, cz-L_outer/2, mesh_size_coarse};
Point(2) = {cx+L_outer/2, cy-L_outer/2, cz-L_outer/2, mesh_size_coarse};
Point(3) = {cx+L_outer/2, cy+L_outer/2, cz-L_outer/2, mesh_size_coarse};
Point(4) = {cx-L_outer/2, cy+L_outer/2, cz-L_outer/2, mesh_size_coarse};
Point(5) = {cx-L_outer/2, cy-L_outer/2, cz+L_outer/2, mesh_size_coarse};
Point(6) = {cx+L_outer/2, cy-L_outer/2, cz+L_outer/2, mesh_size_coarse};
Point(7) = {cx+L_outer/2, cy+L_outer/2, cz+L_outer/2, mesh_size_coarse};
Point(8) = {cx-L_outer/2, cy+L_outer/2, cz+L_outer/2, mesh_size_coarse};

// Rectangle corners
Point(9) = {cx-rect_width/2, cy-rect_height/2, cz, mesh_size_medium};  // Bottom left
Point(10) = {cx+rect_width/2, cy-rect_height/2, cz, mesh_size_medium}; // Bottom right
Point(11) = {cx+rect_width/2, cy+rect_height/2, cz, mesh_size_medium}; // Top right
Point(12) = {cx-rect_width/2, cy+rect_height/2, cz, mesh_size_medium}; // Top left

// Left hole center and points
Point(13) = {cx-d/2, cy, cz, mesh_size_fine}; // Left hole center
Point(14) = {cx-d/2+r, cy, cz, mesh_size_fine}; // Right point of left hole
Point(15) = {cx-d/2, cy+r, cz, mesh_size_fine}; // Top point of left hole
Point(16) = {cx-d/2-r, cy, cz, mesh_size_fine}; // Left point of left hole
Point(17) = {cx-d/2, cy-r, cz, mesh_size_fine}; // Bottom point of left hole

// Right hole center and points
Point(18) = {cx+d/2, cy, cz, mesh_size_fine}; // Right hole center
Point(19) = {cx+d/2+r, cy, cz, mesh_size_fine}; // Right point of right hole
Point(20) = {cx+d/2, cy+r, cz, mesh_size_fine}; // Top point of right hole
Point(21) = {cx+d/2-r, cy, cz, mesh_size_fine}; // Left point of right hole
Point(22) = {cx+d/2, cy-r, cz, mesh_size_fine}; // Bottom point of right hole

// Define the lines (edges of the outer cube)
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// Rectangle edges
Line(13) = {9, 10};  // Bottom edge
Line(14) = {10, 11}; // Right edge
Line(15) = {11, 12}; // Top edge
Line(16) = {12, 9};  // Left edge

// Left hole circles
Circle(17) = {14, 13, 15}; // Left hole quarter 1
Circle(18) = {15, 13, 16}; // Left hole quarter 2
Circle(19) = {16, 13, 17}; // Left hole quarter 3
Circle(20) = {17, 13, 14}; // Left hole quarter 4

// Right hole circles
Circle(21) = {19, 18, 20}; // Right hole quarter 1
Circle(22) = {20, 18, 21}; // Right hole quarter 2
Circle(23) = {21, 18, 22}; // Right hole quarter 3
Circle(24) = {22, 18, 19}; // Right hole quarter 4

// Create curve loops
Curve Loop(1) = {13, 14, 15, 16}; // Rectangle boundary
Curve Loop(2) = {17, 18, 19, 20}; // Left hole
Curve Loop(3) = {21, 22, 23, 24}; // Right hole

// Create rectangle surface with two holes
Plane Surface(1) = {1, 2, 3}; // Rectangle minus both holes
Plane Surface(2) = {2}; // Left hole surface
Plane Surface(3) = {3}; // Right hole surface

// Define the surfaces (faces of the outer cube)
Line Loop(4) = {1, 2, 3, 4};
Plane Surface(4) = {4};

Line Loop(5) = {5, 6, 7, 8};
Plane Surface(5) = {5};

Line Loop(6) = {9, 5, -10, -1};
Plane Surface(6) = {6};

Line Loop(7) = {10, 6, -11, -2};
Plane Surface(7) = {7};

Line Loop(8) = {11, 7, -12, -3};
Plane Surface(8) = {8};

Line Loop(9) = {12, 8, -9, -4};
Plane Surface(9) = {9};

// Create cube volume with embedded rectangle surfaces
Surface Loop(1) = {4, 5, 6, 7, 8, 9, 1, 2, 3}; // Include rectangle surfaces
Volume(1) = {1};

// Physical groups
Physical Volume("domain", 1) = {1};
Physical Surface("cube_wall_m_xy", 2) = {4}; // Outer bottom
Physical Surface("cube_wall_p_xy", 3) = {5}; // Outer top
Physical Surface("cube_wall_m_xz", 4) = {6}; // Outer front
Physical Surface("cube_wall_p_yz", 5) = {7}; // Outer right
Physical Surface("cube_wall_p_xz", 6) = {8}; // Outer back
Physical Surface("cube_Wall_m_yz", 7) = {9}; // Outer left
Physical Surface("rectangle_surface", 8) = {1};
Physical Surface("left_hole_surface", 9) = {2};
Physical Surface("right_hole_surface", 10) = {3};

// Mesh size fields for better control
Field[1] = Distance;
Field[1].CurvesList = {17, 18, 19, 20}; // Left hole edges
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = mesh_size_fine;
Field[2].SizeMax = mesh_size_medium;
Field[2].DistMin = 0.2;
Field[2].DistMax = 1.0;

Field[3] = Distance;
Field[3].CurvesList = {21, 22, 23, 24}; // Right hole edges
Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = mesh_size_fine;
Field[4].SizeMax = mesh_size_medium;
Field[4].DistMin = 0.2;
Field[4].DistMax = 1.0;

Field[5] = Distance;
Field[5].CurvesList = {13, 14, 15, 16}; // Rectangle edges
Field[6] = Threshold;
Field[6].InField = 5;
Field[6].SizeMin = mesh_size_medium;
Field[6].SizeMax = mesh_size_coarse;
Field[6].DistMin = 0.5;
Field[6].DistMax = 2.5;

// Field to make cube surfaces coarse
Field[7] = Constant;
Field[7].VIn = mesh_size_coarse;
Field[7].VOut = mesh_size_coarse;
Field[7].SurfacesList = {4, 5, 6, 7, 8, 9};

Field[8] = Min;
Field[8].FieldsList = {2, 4, 6};
Field[9] = Max;
Field[9].FieldsList = {8, 7};
Background Field = 9;

// Mesh optimization
Mesh.OptimizeNetgen = 1;
Mesh.Optimize = 1;
Mesh.ElementOrder = 1;