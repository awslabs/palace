// Define a unit cube with 3 elements along each direction

// Define the size of the cube and number of elements
lc = 1.0 / 3.0; // characteristic length (edge length of the smallest elements)

// Define points
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};
Point(5) = {0, 0, 1, lc};
Point(6) = {1, 0, 1, lc};
Point(7) = {1, 1, 1, lc};
Point(8) = {0, 1, 1, lc};

// Define lines
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

// Define surfaces
Line Loop(1) = {1, 2, 3, 4}; // bottom face
Plane Surface(1) = {1};

Line Loop(2) = {5, 6, 7, 8}; // top face
Plane Surface(2) = {2};

Line Loop(3) = {1, 10, -5, -9}; // front face
Plane Surface(3) = {3};

Line Loop(4) = {2, 11, -6, -10}; // back face
Plane Surface(4) = {4};

Line Loop(5) = {3, 12, -7, -11}; // right face
Plane Surface(5) = {5};

Line Loop(6) = {4, 9, -8, -12}; // left face
Plane Surface(6) = {6};

// Define volume
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// Define physical entities
Physical Volume(1) = {1}; // volume
Physical Surface(1) = {1}; // bottom face
Physical Surface(2) = {2}; // top face
Physical Surface(3) = {3}; // front face
Physical Surface(4) = {4}; // back face
Physical Surface(5) = {5}; // right face
Physical Surface(6) = {6}; // left face

// Meshing options
Transfinite Line {1, 2, 3, 4, 5, 6, 7, 8} = 3; // divide each edge into 3 segments

// Generate mesh
Mesh.RecombineAll = 1; // ensure elements are tetrahedra (if using 3D)

// Save mesh
Mesh 3;
