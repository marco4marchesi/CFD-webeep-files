
// Flat plat mesh


X1 = -0.06;
X2 = 0.0;
X3 = 0.3048;
Y1 = 0;
Y2 = 0.03;


//Points
Point(1)={X1, Y1, 0};
Point(2)={X2, Y1, 0};
Point(3)={X3, Y1, 0};
Point(4)={X3, Y2, 0};
Point(5)={X2, Y2, 0};
Point(6)={X1, Y2, 0};

// Lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};

// Construction line
Line(7) = {5,2};


//Transfinite curve
Transfinite Curve{1}= 60 Using Progression 0.95;
Transfinite Curve{5}= 60 Using Progression 1/0.95;

Transfinite Curve{2}= 100 Using Progression 1.05;
Transfinite Curve{4}= 100 Using Progression 1/1.05;

Transfinite Curve{6}= 64 Using Progression 0.9;
Transfinite Curve{7}= 64 Using Progression 0.9;
Transfinite Curve{3}= 64 Using Progression 1/0.9;


// Build surfaces
Curve Loop(11)={1,-7,5,6};
Plane Surface(1) = {11};
Transfinite Surface{1}={1,2,5,6};
Recombine Surface{1};

Curve Loop(12)={2,3,4,7};
Plane Surface(2) = {12};
Transfinite Surface{2}={2,3,4,5};
Recombine Surface{2};


//Physical names
Physical Surface("VOLUME")={1,2};
Physical Line("symmetry")={1};
Physical Line("wall")={2};
Physical Line("outlet")={3};
Physical Line("farfield")={4,5};
Physical Line("inlet")={6};





