Mesh.Algorithm=5;
Point(1)={934000,20000,0.0,5000};
Point(2)={934000,-20000,0.0,5000};
Point(3)={1174000,-20000,0.0,5000};
Point(4)={1174000,20000,0.0,5000};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
Line Loop(5)={1,2,3,4};
Physical Line(6)={1};
Physical Line(7)={2};
Physical Line(8)={3};
Physical Line(9)={4};
Plane Surface(10)={5};
Physical Surface(11)={10};
Transfinite Surface {10};
Recombine Surface {10};