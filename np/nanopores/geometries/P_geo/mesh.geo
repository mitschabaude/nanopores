/*****
gmsh geo for hourglass-shaped pore like in PRE
*****/

nm = 1e-9;
// pore radius min @center, max
r0 = 5*nm;
r1 = 20*nm;
// pore length min @center, membrane thickness
l0 = 20*nm;
l1 = 50*nm;
// Radius of Omega
Ry = 100*nm;
Rx = 25*nm;
// DNA size
DNAlength = 100*nm; // one Kuhn-length
DNAradius = 1.1*nm;

// characteristic length
lc = nm*5;
lcDNA = lc/2; //*6e-1;
lcFluid = lc ;
lcMembrane = lc/2; //*8e-1;

// Fluid
Point(1) = {0, Ry, 0, lcFluid};
Point(2) = {Rx, Ry, 0, lcFluid};
Point(3) = {Rx, -Ry, 0, lcFluid};
Point(4) = {0, -Ry, 0, lcFluid};

// DNA
Point(5) = {0, DNAlength/2.0, 0, lcMembrane};
Point(6) = {DNAradius, DNAlength/2.0, 0, lcDNA};
Point(7) = {DNAradius, -DNAlength/2.0, 0, lcDNA};
Point(8) = {0, -DNAlength/2.0, 0, lcMembrane};

// Membrane
Point(9) = {r0, -l0/2.0, 0, lcMembrane};
Point(10) = {r0, l0/2.0, 0, lcMembrane};
Point(11) = {r1, l1/2.0, 0, lcMembrane};
Point(12) = {Rx, l1/2.0, 0, lcMembrane};
Point(13) = {Rx, -l1/2.0, 0, lcMembrane};
Point(14) = {r1, -l1/2.0, 0, lcMembrane};

// Line
Line(1) = {5, 1};
Line(2) = {1, 2};
Line(3) = {2, 12};
Line(4) = {12, 11};
Line(5) = {11, 10};
Line(6) = {10, 9};
Line(7) = {9, 14};
Line(8) = {14, 13};
Line(9) = {13, 3};
Line(10) = {3, 4};
Line(11) = {4, 8};
Line(12) = {8, 7};
Line(13) = {7, 6};
Line(14) = {6, 5};

// Fluid
Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
Plane Surface(1) = {1};

// DNA
Line(15) = {8, 5};
Line Loop(2) = {-12, -13, -14, 15};
Plane Surface(2) = {2};

// Membrane
Line(16) = {12, 13};
Line Loop(3) = {-4, -5, -6, -7, -8, 16};
Plane Surface(3) = {3};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};

// We could also use a Box field to impose a step change in element
// sizes inside a box
Field[4] = Box;
Field[4].VIn = lcDNA*0.75;
Field[4].VOut = lc;
Field[4].XMin = DNAradius*0.75;
Field[4].XMax = DNAradius*1.5;
Field[4].YMin = -DNAlength*0.5 - DNAradius*0.5;
Field[4].YMax = DNAlength*0.5 + DNAradius*0.5;

// We could also use a Box field to impose a step change in element
// sizes inside a box
Field[5] = Box;
Field[5].VIn = lcDNA ;
Field[5].VOut = lc;
Field[5].XMin = r0 - DNAradius*0.5 ;
Field[5].XMax = r0 + DNAradius*0.25 ;
Field[5].YMin = -l0*0.5 ;
Field[5].YMax = l0*0.5 ;

Field[7] = Min;
Field[7].FieldsList = {4, 5};

// Use background fields ?
//Background Field = 7;

// Don't extend the elements sizes from the boundary inside the domain
Mesh.CharacteristicLengthExtendFromBoundary = 1;

// 2D mesh algorithm (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal,7=bamg, 8=delquad)
Mesh.Algorithm = 5;

// Disable question dialogs
General.ExpertMode = 1;
