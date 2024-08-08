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

Periodic Surface {out[4]} = {out[2]} Translate {0, 1, 0};
Periodic Surface {out[3]} = {out[5]} Translate {1, 0, 0};

Physical Surface(1) = {out[3]}; // left (x-direction)
Physical Surface(2) = {out[5]}; // right (x-direction)
Physical Surface(3) = {out[2]}; // bottom (y-direction)
Physical Surface(4) = {out[4]}; // top (y-direction)
Physical Surface(5) = {1};      // back (z-direction)
Physical Surface(6) = {out[0]}; // front (z-direction)
