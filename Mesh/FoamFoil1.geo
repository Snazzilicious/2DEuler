Mesh.SaveAll=1;

//Domain Boundary Points
//Bottom left
Point(1) = {-4.0, -3.0, 0, 1.0};
//Top left
Point(2) = {-4.0, 3.0, 0, 1.0};
//Top right
Point(3) = {10.0, 3.0, 0, 1.0}; 
//Bottom right
Point(4) = {10.0, -3.0, 0, 1.0}; 


// Domain boundary lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};




// Airfoil point definitions
thickness = 0.1;
overlap1 = 0.4;
overlap2 = 0.4;
len1 = 1.0;
len2 = 2.0;
len3 = 1.75;

scale = 0.03;

// Bottom
Point(5) = {0.0, 0.0, 0.0, scale};
Point(6) = {len1, 0.0, 0.0, scale};
Point(7) = {len1, thickness, 0.0, scale};
Point(8) = {len1+len2-overlap1-overlap2, thickness, 0.0, scale};
Point(9) = {len1+len2-overlap1-overlap2, 0.0, 0.0, scale};
Point(10) = {len1+len2-overlap1-overlap2+len3, 0.0, 0.0, scale};
//top
Point(11) = {0.0, thickness, 0.0, scale};
Point(12) = {len1 - overlap1, thickness, 0.0, scale};
Point(13) = {len1 - overlap1, 2*thickness, 0.0, scale};
Point(14) = {len1 - overlap1 + len2, 2*thickness, 0.0, scale};
Point(15) = {len1 - overlap1 + len2, thickness, 0.0, scale};
Point(16) = {len1 - overlap1 + len2 + len3 - overlap2, thickness, 0.0, scale};


////body
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};

Line(10) = {10,16};

Line(11) = {16,15};
Line(12) = {15,14};
Line(13) = {14,13};
Line(14) = {13,12};
Line(15) = {12,11};

Line(16) = {11,5};

// Identify domain boundaries
Curve Loop(1) = {1,2,3,4, 5,6,7,8,9,10,11,12,13,14,15,16};
Plane Surface(1) = {1};

//// Body
// Upper and lower boundary
Transfinite Line {-2,4} = 27 Using Progression .97;
// forward and backward boundary
Transfinite Line {1} = 19;
Transfinite Line {-3} = 11;


Physical Curve("Inflow") = {1};
Physical Curve("Outflow") = {3};
Physical Curve("Top") = {2};
Physical Curve("Bottom") = {4};
Physical Curve("Body") = {5,6,7,8,9,10,11,12,13,14,15,16};




