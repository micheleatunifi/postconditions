% Matlab script for safety of semialgebraic continuous system via SOS programming. See 
% "Complete algorithms for algebraic strongest postconditions and
% weakest preconditions in  polynomial ODEs" (Section 7)
% by Michele Boreale, 2019.
% Requires SOSTOOLS and SeDuMi (or an equivalent tool for semidefinite programming)


% 3D LOTKA-VOLTERRA

clear;
pvar x0 y0 x y z;
PLV = sosprogram([x0; y0; x; y; z]);
mon2 = monomials([x0; y0; x; y; z], [1]);
mon1 = monomials([x0; y0; x; y; z],[0 1]);

[PLV, s1]=sossosvar(PLV,mon2);
[PLV, s2]=sossosvar(PLV,mon2);
[PLV, t1]=sospolyvar(PLV,mon1);
[PLV, t2]=sospolyvar(PLV,mon1);

pin = -((x0-2)^2+(y0-2)^2)+1.15^2;
pu  =-((x-.5)^2+(y-5)^2)+1.5^2;
p1=x*z*y-x*3*y0-z*3*y0-y*3*y0+3^2*y0+3*y0^2;
p2=x0+y0+3-x-y-z;
c=-(s1*pin+s2*pu+t1*p1+t2*p2+1);

PLV=sosineq(PLV,c);
PLV=sossolve(PLV);

S1=sosgetsol(PLV,s1);
S2=sosgetsol(PLV,s2);
T1=sosgetsol(PLV,t1);
T2=sosgetsol(PLV,t2);

CLV=-(S1*pin+S2*pu+T1*p1+T2*p2+1); 


% COUPLED SPRING MASS SYSTEM

clear;
pvar x10 x20 x1 x2 v1 v2;
PSM = sosprogram([x10; x20; x1; x2; v1; v2 ]);
mon2 = monomials([x10; x20; x1; x2; v1; v2 ], [0 1]);
mon1 = monomials([x10; x20; x1; x2; v1; v2 ], [ 1]);

[PSM, s1]=sossosvar(PSM,mon2);
[PSM, s2]=sossosvar(PSM,mon2);
[PSM, t1]=sospolyvar(PSM,mon1);
[PSM, t2]=sospolyvar(PSM,mon1);

pin1 = -(x10-1/2)^2-(x20-3/2)^2+1/4;
pu = x2- x1-2.17;
q1 = -v1^2 + 2*v1*v2 - 3*x1^2 + 4*x1*x2 + 3*x10^2 - 4*x10*x20 - x2^2 + x20^2 - 2*x1 + 2*x10;
q2 = v1^2 + 4*v1*v2 + 3*v2^2 + 2*x1*x2 - 2*x10*x20 + x2^2 - x20^2 - 4*x1 + 4*x10 - 6*x2 + 6*x20;
c= -(s1*pin1 + s2*pu  + t1*q1+ t2*q2 + 1); 

PSM=sosineq(PSM,c);
PSM=sossolve(PSM);

S1=sosgetsol(PSM,s1);
S2=sosgetsol(PSM,s2);
T1=sosgetsol(PSM,t1);
T2=sosgetsol(PSM,t2);

CSM= -(S1*pin1 +   S2*pu   + T1*q1+ T2*q2+1);


