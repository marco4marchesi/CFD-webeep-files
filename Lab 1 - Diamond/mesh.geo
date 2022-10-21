// DIAMOND AIRFOIL MESH

R = 30;
c = 1;
t = 0.3;
h = 0.01;
H = 2;


//========================= POINTS

//airfoil
Point(1) = {0, 0, 0, h/5};
Point(2) = {c/2, t/2, 0, h};
Point(3) = {c, 0, 0, h/5};
Point(4) = {c/2, -t/2, 0, h};

//farfield
Point(5) = {R+15, 0, 0, H};
Point(6) = {15, R, 0, H};
Point(7) = {-R+15, 0, 0, H};
Point(8) = {15, -R, 0, H};

// farfield center
Point(9) = {15,0,0,h};

//========================== LINES

// airfoil
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

// farfield
Circle(5) = {5, 9, 6};
Circle(6) = {6, 9, 7};
Circle(7) = {7, 9, 8};
Circle(8) = {8, 9, 5};



//=========================== LOOPS
//airfoil
Line Loop(1) = {1,2,3,4};
//farfield
Line Loop(2) = {5,6,7,8};

//=========================== SURFACES
Plane Surface(1) = {1,2};

//===============================MESH
Physical Surface("Volume") = {1};
Physical Line("AIRFOIL") = {1,2,3,4};
Physical Line("FARFIELD") = {5,6,7,8};

Mesh.RandomFactor = 1e-12;
//Mesh.Algorithm = 1;
