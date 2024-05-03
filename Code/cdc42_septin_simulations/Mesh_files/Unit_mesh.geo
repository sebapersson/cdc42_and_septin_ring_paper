// The parameter "dens_fine" which determines the number
// of nodes
dens_fine = 0.045; // Corresponding to 8964 nodes

// Macro which creates a sphere with provided center 
// coordinates, radious and surface labels. 
// Args:
//    x/y/z_center, center coordinates 
//    rad, the radious of the sphere 
//    i_surf, i_phys, the indices for the surface 
Macro MakeSphere

p1 = newp; Point(p1) = {x_cent, y_cent, z_cent, ndens};
p2 = newp; Point(p2) = {x_cent + rad, y_cent, z_cent, ndens};
p3 = newp; Point(p3) = {x_cent-rad, y_cent, z_cent, ndens};
p4 = newp; Point(p4) = {x_cent, y_cent+rad, z_cent, ndens};
p5 = newp; Point(p5) = {x_cent, y_cent-rad, z_cent, ndens};
p6 = newp; Point(p6) = {x_cent, y_cent, z_cent+rad, ndens};
p7 = newp; Point(p7) = {x_cent, y_cent, z_cent-rad, ndens};

c1 = newc; Circle(c1) = {p7, p1, p4};
c2 = newc; Circle(c2) = {p2, p1, p4};
c3 = newc; Circle(c3) = {p6, p1, p4};
c4 = newc; Circle(c4) = {p3, p1, p4};
c5 = newc; Circle(c5) = {p5, p1, p2};
c6 = newc; Circle(c6) = {p7, p1, p5};
c7 = newc; Circle(c7) = {p5, p1, p6};
c8 = newc; Circle(c8) = {p3, p1, p5};

l1 = newll; Curve Loop(l1) = {c1, -c2, -c5, -c6};
l2 = newll; Curve Loop(l2) = {c3, -c2, -c5, c7};
l3 = newll; Curve Loop(l3) = {c4, -c3, -c7, -c8};
l4 = newll; Curve Loop(l4) = {c1, -c4, c8, -c6};

s1 = news; Surface(s1) = {l1};
s2 = news; Surface(s2) = {l2};
s3 = news; Surface(s3) = {l3};
s4 = news; Surface(s4) = {l4};


Surface Loop(i_surf) = {s1, s4, s3, s2};
Physical Surface(i_phys) = {s1, s3, s2, s4};

Return

ndens = 3e-7
x_cent = 0.0; y_cent = 0.0; z_cent = 0.0;
rad = 1.0;
i_surf = 1; i_phys = 1;
Call MakeSphere;

// The actual volume 
Volume(1) = {1};
Physical Volume(1) = {1};

// Define Ball field
Field[1] = Ball;
Field[1].Radius = 0.9;
Field[1].Thickness = 0.0; //Only for later versions of gmsh
Field[1].VIn = 0.4;
Field[1].VOut = dens_fine;

Background Field = 1;

// Do not use the boundary cell size in the interior
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;


