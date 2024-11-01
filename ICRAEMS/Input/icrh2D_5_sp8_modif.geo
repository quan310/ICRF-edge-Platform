xmin= -1.;
xmax= 1.;
ymin= 0.;
ymax= 1.;
nx1=100;
nx2=3*nx1;
ny=nx1;
lc1 = (xmax-xmin)/nx1;
lc2 = (xmax-xmin)/nx2;
hpml= 0.15;
npml= 9;
Sc=3;
mxpml1 = nx1;
mxpml2 = nx1*Sc;
mypml = ny;
ra=1+((ny*ny-2*ny*(ny-1)*(1-Sc))^0.5-ny)/(ny*(ny-1));

antwidth= 0.2;
antthickness= 0.02;
antdx= 0.14;
anthigh= 0.18;

centx=(xmax-xmin)/2+xmin;
antx= antdx/2;
antx2= antx+antwidth;
anty= anthigh+ymin;
anty2= anty+antthickness;
Point(1) = {xmin,ymin,0,lc1};
Point(2) = {xmin+hpml,ymin,0,lc1};
Point(3) = {xmax-hpml,ymin,0,lc1};
Point(4) = {xmax,ymin,0,lc1};
Point(5) = {xmin,ymin+hpml,0,lc1};
Point(6) = {xmin+hpml,ymin+hpml,0,lc1};
Point(7) = {xmax-hpml,ymin+hpml,0,lc1};
Point(8) = {xmax,ymin+hpml,0,lc1};

Point(9) = {xmin,ymax-hpml,0,lc1};
Point(10) = {xmin+hpml,ymax-hpml,0,lc1};
Point(11) = {xmax-hpml,ymax-hpml,0,lc1};
Point(12) = {xmax,ymax-hpml,0,lc1};
Point(13) = {xmin,ymax,0,lc1};
Point(14) = {xmin+hpml,ymax,0,lc1};
Point(15) = {xmax-hpml,ymax,0,lc1};
Point(16) = {xmax,ymax,0,lc1};

Point(17) = {centx-antx2,anty,0,lc2};
Point(18) = {centx-antx,anty,0,lc2};
Point(19) = {centx-antx,anty2,0,lc2};
Point(20) = {centx-antx2,anty2,0,lc2};
Point(21) = {centx+antx,anty,0,lc2};
Point(22) = {centx+antx2,anty,0,lc2};
Point(23) = {centx+antx2,anty2,0,lc2};
Point(24) = {centx+antx,anty2,0,lc2};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {5,6};
Line(5) = {6,7};
Line(6) = {7,8};
Line(7) = {9,10};
Line(8) = {10,11};
Line(9) = {11,12};
Line(10) = {13,14};
Line(11) = {14,15};
Line(12) = {15,16};
Line(13) = {1,5};
Line(14) = {2,6};
Line(15) = {3,7};
Line(16) = {4,8};
Line(17) = {5,9};
Line(18) = {6,10};
Line(19) = {7,11};
Line(20) = {8,12};
Line(21) = {9,13};
Line(22) = {10,14};
Line(23) = {11,15};
Line(24) = {12,16};

Line(25) = {17,18};
Line(26) = {18,19};
Line(27) = {19,20};
Line(28) = {20,17};
Line(29) = {21,22};
Line(30) = {22,23};
Line(31) = {23,24};
Line(32) = {24,21};

Transfinite Curve {13,14,15,16,21,22,23,24} = npml Using Progression 1;
Transfinite Curve {1,4,7,10,3,6,9,12} = npml Using Progression 1;
Transfinite Curve {8,11} = mxpml1 Using Progression 1;
Transfinite Curve {2,5}  = mxpml2 Using Progression 1;
Transfinite Curve {17,18,19,20} = mypml Using Progression ra;

Curve Loop(1) = {1,14,-4,-13};
Plane Surface(1) = {1};
Transfinite Surface{1};
Recombine Surface{1};

Curve Loop(2) = {2,15,-5,-14};
Plane Surface(2) = {2};
Transfinite Surface{2};
Recombine Surface{2};

Curve Loop(3) = {3,16,-6,-15};
Plane Surface(3) = {3};
Transfinite Surface{3};
Recombine Surface{3};

Curve Loop(4) = {4,18,-7,-17};
Plane Surface(4) = {4};
Transfinite Surface{4};
Recombine Surface{4};

Curve Loop(5) = {5,19,-8,-18};
Curve Loop(10)= {25,26,27,28};
Curve Loop(11)= {29,30,31,32};
Plane Surface(5) = {5,10,11};
Plane Surface(10)= {10};
Plane Surface(11)= {11};

Curve Loop(6) = {6,20,-9,-19};
Plane Surface(6) = {6};
Transfinite Surface{6};
Recombine Surface{6};

Curve Loop(7) = {7,22,-10,-21};
Plane Surface(7) = {7};
Transfinite Surface{7};
Recombine Surface{7};

Curve Loop(8) = {8,23,-11,-22};
Plane Surface(8) = {8};
Transfinite Surface{8};
Recombine Surface{8};

Curve Loop(9) = {9,24,-12,-23};
Plane Surface(9) = {9};
Transfinite Surface{9};
Recombine Surface{9};


