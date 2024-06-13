// Periodic hexahedral cube mesh

Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};

Characteristic Length {:} = 0.25;

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//Periodic Curve {1} = {-3};
//Periodic Curve {2} = {-4};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface {1};

Recombine Surface {1};
out[] = Extrude {0, 0, 1} { Surface{1}; Layers{4}; Recombine; };

Physical Volume(1) = {out[1]};

Mesh 3;
Mesh.MshFileVersion = 2.2;

//Periodic Surface {out[0]} = {1} Translate {0, 0, 1}; // No periodicity in the z-direction (back & front)
Periodic Surface {out[4]} = {out[2]} Translate {0, 1, 0}; // left and right
Periodic Surface {out[3]} = {out[5]} Translate {1, 0, 0}; // bottom and top

Physical Surface(1) = {1};      // back
Physical Surface(2) = {out[5]}; // front
Physical Surface(3) = {out[4]}; // left
Physical Surface(4) = {out[2]}; // right
Physical Surface(5) = {out[0]}; // bottom
Physical Surface(6) = {out[3]}; // top
