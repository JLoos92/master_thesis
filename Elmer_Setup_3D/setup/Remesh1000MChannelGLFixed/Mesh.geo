Mesh.Algorithm=5;
Point(1)={1010000,5000,0.0,750};
Point(2)={1010000,-5000,0.0,750};
Point(3)={1080000,-5000,0.0,750};
Point(4)={1080000,5000,0.0,750};
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
Field[1] = Box;
Field[1].VIn = 100;
Field[1].VOut = 750;
Field[1].XMax = 1080000;
Field[1].XMin = 1053500;
Field[1].YMax = 700;
Field[1].YMin = -700;
Field[1].ZMax = 0;
Field[1].ZMin = 0;
Background Field = 1;
