// Define the cube sizes
L_outer = 15.0;
R = 3.0;
r = 1.0;

// Set mesh sizes for different regions
mesh_size_coarse = 2.0;  // Cube surfaces (very coarse)
mesh_size_medium = 0.4;  // Disk region
mesh_size_fine = 0.2;    // Hole region

// Mesh algorithm and control
Mesh.Algorithm3D = 1; // Delaunay algorithm for 3D mesh
Mesh.CharacteristicLengthFactor = 1.0;
Mesh.MshFileVersion = 2.2;
Mesh.CharacteristicLengthMin = 0.05;
Mesh.CharacteristicLengthMax = 3.0;

// Define center point for the cube and the disk
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

Point(9) = {cx, cy, cz, mesh_size_medium}; // Center point
Point(10) = {cx+R, cy, cz, mesh_size_medium}; // Point on outer circle 
Point(11) = {cx, cy+R, cz, mesh_size_medium}; // Point on outer circle
Point(12) = {cx-R, cy, cz, mesh_size_medium}; // Point on outer circle
Point(13) = {cx, cy-R, cz, mesh_size_medium}; // Point on outer circle
Point(14) = {cx+r, cy, cz, mesh_size_fine}; // Point on inner circle
Point(15) = {cx, cy+r, cz, mesh_size_fine}; // Point on inner circle
Point(16) = {cx-r, cy, cz, mesh_size_fine}; // Point on inner circle
Point(17) = {cx, cy-r, cz, mesh_size_fine}; // Point on inner circle

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


Circle(13) = {10, 9, 11}; // Outer circle quarter 1
Circle(14) = {11, 9, 12}; // Outer circle quarter 2
Circle(15) = {12, 9, 13}; // Outer circle quarter 3
Circle(16) = {13, 9, 10}; // Outer circle quarter 4
Circle(17) = {14, 9, 15}; // Inner circle quarter 1
Circle(18) = {15, 9, 16}; // Inner circle quarter 2
Circle(19) = {16, 9, 17}; // Inner circle quarter 3
Circle(20) = {17, 9, 14}; // Inner circle quarter 4

// Create curve loops
Curve Loop(1) = {13, 14, 15, 16}; // Outer circle
Curve Loop(2) = {17, 18, 19, 20}; // Inner circle (hole)

// Create disk surface with hole
Plane Surface(1) = {1, 2}; // Outer circle minus inner circle
Plane Surface(2) = {2}; // Hole surface

// Create cube surfaces
//Curve Loop(3) = {1, 2, 3, 4}; // Bottom face
//Curve Loop(4) = {5, 6, 7, 8}; // Top face
//Curve Loop(5) = {1, 10, -5, -9}; // Front face
//Curve Loop(6) = {3, 12, -7, -11}; // Back face
//Curve Loop(7) = {2, 11, -6, -10}; // Right face
//Curve Loop(8) = {4, 9, -8, -12}; // Left face

// Define the surfaces (faces of the outer cube)
Line Loop(3) = {1, 2, 3, 4};
Plane Surface(3) = {3};

Line Loop(4) = {5, 6, 7, 8};
Plane Surface(4) = {4};

Line Loop(5) = {9, 5, -10, -1};
Plane Surface(5) = {5};

Line Loop(6) = {10, 6, -11, -2};
Plane Surface(6) = {6};

Line Loop(7) = {11, 7, -12, -3};
Plane Surface(7) = {7};

Line Loop(8) = {12, 8, -9, -4};
Plane Surface(8) = {8};

//Plane Surface(2) = {3}; // Bottom
//Plane Surface(3) = {4}; // Top
//Plane Surface(4) = {5}; // Front
//Plane Surface(5) = {6}; // Back
//Plane Surface(6) = {7}; // Right
//Plane Surface(7) = {8}; // Left

// Create cube volume with embedded disk surfaces
Surface Loop(1) = {3, 4, 5, 6, 7, 8, 1, 2}; // Include disk surfaces
Volume(1) = {1};

// Physical groups
Physical Volume("domain", 1) = {1};
Physical Surface("cube_wall_m_xy", 2) = {3}; // Outer bottom
Physical Surface("cube_wall_p_xy", 3) = {4}; // Outer top
Physical Surface("cube_wall_m_xz", 4) = {5}; // Outer front
Physical Surface("cube_wall_p_yz", 5) = {6}; // Outer right
Physical Surface("cube_wall_p_xz", 6) = {7}; // Outer back
Physical Surface("cube_Wall_m_yz", 7) = {8}; // Outer left
Physical Surface("disk_surface", 8) = {1};
Physical Surface("hole_surface", 9) = {2};
//Physical Curve("disk_edge", 5) = {13, 14, 15, 16};
//Physical Curve("hole_edge", 6) = {17, 18, 19, 20};


// Mesh size fields for better control
Field[1] = Distance;
Field[1].CurvesList = {17, 18, 19, 20}; // Hole edges
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = mesh_size_fine;
Field[2].SizeMax = mesh_size_medium;
Field[2].DistMin = 0.2;
Field[2].DistMax = 1.0;

Field[3] = Distance;
Field[3].CurvesList = {13, 14, 15, 16}; // Disk edges
Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = mesh_size_medium;
Field[4].SizeMax = mesh_size_coarse;
Field[4].DistMin = 0.5;
Field[4].DistMax = 2.5;

// Field to make cube surfaces coarse
Field[6] = Constant;
Field[6].VIn = mesh_size_coarse;
Field[6].VOut = mesh_size_coarse;
Field[6].SurfacesList = {3, 4, 5, 6, 7, 8};

Field[5] = Min;
Field[5].FieldsList = {2, 4};
Field[7] = Max;
Field[7].FieldsList = {5, 6};
Background Field = 7;

// Mesh optimization
Mesh.OptimizeNetgen = 1;
Mesh.Optimize = 1;
Mesh.ElementOrder = 1;