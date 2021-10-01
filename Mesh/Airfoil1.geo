Mesh.SaveAll=1;

//Domain Boundary Points
//Bottom left
Point(1) = {-4.0, -3.0, 0, 1.0};
//Top left
Point(2) = {-4.0, 3.0, 0, 1.0};
//Top right
Point(3) = {6.0, 3.0, 0, 1.0}; 
//Bottom right
Point(4) = {6.0, -3.0, 0, 1.0}; 


// Domain boundary lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};



//Airfoil points
Point(5) = {-1.0, 0.0, 0, 1.0}; //leading edge
Point(6) = {2.0, 0.0, 0, 1.0}; // trailing edge
// center point
Point(7) = {0.0, 0.0 , 0 , 1.0};

// crossover pts
Point(8) = {0.5, 0.2165 , 0 , 1.0}; //upper
Point(9) = {0.5, -0.1299 , 0 , 1.0}; //lower



//body
Ellipse(5) = {5, 7, 5, 8};
Ellipse(6) = {5, 7, 5, 9};
Line(7) = {8,6};
Line(8) = {9,6};

// Identify domain boundaries
Curve Loop(1) = {1,2,3,4,5,7,-8,-6};
Plane Surface(1) = {1};

// Body
Transfinite Line {5} = 11;
Transfinite Line {6} = 11;
Transfinite Line {7} = 11;
Transfinite Line {8} = 11;
// Upper and lower boundary
Transfinite Line {-2,4} = 11 Using Progression .97;
// forward and backward boundary
Transfinite Line {1} = 11;
Transfinite Line {-3} = 6;


Physical Curve("Inflow") = {1};
Physical Curve("Outflow") = {3};
Physical Curve("Top") = {2};
Physical Curve("Bottom") = {4};
Physical Curve("Body") = {5,6,7,8};




