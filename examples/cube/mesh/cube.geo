SetFactory("OpenCASCADE");
gridsize = 100/3;

// Face
// -- Points
Point(1) = {0,0,0,gridsize};
Point(2) = {0,0,100,gridsize};
Point(3) = {0,100,100,gridsize};
Point(4) = {0,100,0,gridsize};

// -- Line
Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};
Line Loop(9) = {5,6,7,8};

// -- Surface
Plane Surface(10) = {9};
Recombine Surface{10};

// Extrude to 3D
bodyExtrusion[] =
Extrude { 100,0,0 }
{
    Surface{10};
    Layers{100/gridsize};
    Recombine;
};

// Boundaries
Physical Surface("back") = {10};
Physical Surface("right") = {bodyExtrusion[2]};
Physical Surface("top") = {bodyExtrusion[3]};
Physical Surface("left") = {bodyExtrusion[4]};
Physical Surface("bottom") = {bodyExtrusion[5]};
Physical Surface("front") = {bodyExtrusion[0]};

// Volume
Physical Volume("mesh") = {bodyExtrusion[1]};
Mesh 3;
Coherence Mesh;