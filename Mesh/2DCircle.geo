Mesh.SaveAll=1;

//Domain Boundary Points
//Bottom left
Point(1) = {0, 0, 0, 1.0};
//Top left
Point(2) = {0, 1.0, 0, 1.0};
//Top right
Point(3) = {2.0, 1.0, 0, 1.0}; 
//Bottom right
Point(4) = {2.0, 0, 0, 1.0}; 


// Domain boundary lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};



//Circle points
radius = .25;
Point(5) = {0.5, 0.6, 0, 1.0};
Point(6) = {0.5, 0.4, 0, 1.0};
//circle center point
Point(7) = {0.5, 0.5 , 0 , 1.0};

//body
Circle(5) = {6, 7, 5};
Circle(6) = {5, 7, 6};


// Identify domain boundaries
Curve Loop(1) = {1,2,3,4,-5,-6};
Plane Surface(1) = {1};

// Body
Transfinite Line {5,6} = 11;
// Upper and lower boundary
Transfinite Line {-2,4} = 11 Using Progression .97;
// forward and backward boundary
Transfinite Line {1} = 11;
Transfinite Line {-3} = 6;


Physical Curve("Inflow") = {1};
Physical Curve("Outflow") = {3};
Physical Curve("Top") = {2};
Physical Curve("Bottom") = {4};
Physical Curve("Body") = {5,6};
