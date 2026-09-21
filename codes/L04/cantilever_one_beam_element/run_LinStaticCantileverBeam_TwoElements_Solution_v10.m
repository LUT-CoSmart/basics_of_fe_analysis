% BK10A6400 Basics of FE Analysis (FEMBasics2025)
% Teacher in charge: Marko Matikainen (MKM)
% Example: Cantilever problem
% Solution for the simple cantilever problem.
%
% Goal: The code solves global displacements at node 3 and member forces of a simple beam structure. 
% Coded by MKM for student's usage in the FEMBasics2025 course

%          N1              N2
%          L1,I1     L2,I2
%          \               |F 
%          \o-----o-------o
%          \  (1)     (2) 
%
%  v1,theta1           v2,theta2 

% So-called connectivity matrix
% Ele1 N1 N2   v1,theta1,v2,theta2


% First matlab-based FE code for a simple beam structure
clear all;
close all;
clc;
%format long;
format shortG;

% Units are mm and N and radians
% Beam structure
E=207000;       % Young's modulus
L=2000;         % length
H=100;          % height 
W=50;           % width
Iz=W*H^3/12;
alpha=0;

% Force (in vertical (global) direction at node 2)
F=-1000;

% Elemental stiffness matrices
Kloc1=KlocEB4DOF(E,Iz,L/2);
Kloc2=KlocEB4DOF(E,Iz,L/2);


% Transformation into a global coordinate system
% Transformation matrix for a beam 1 
T1=TEB4DOFs(alpha);

% Global elemental stiffness matrix for a rod 1
Kglob1=T1'*Kloc1*T1;
Kglob2=T1'*Kloc2*T1;

% As you can see, local stiffness matrices have zero rows and columns due
% to transformation to global coordinate system and paraller coordinate
% systems so no point to make transformation because zero rows and columns 
% need to eliminated. Note that this is just special case. For 6 DOFs beam
% element   later, use transformation.

%Let's use Klocs directly:
Kglob=zeros(6,6);
Kglob(1:4,1:4)=Kloc1;
Kglob(3:6,3:6)=Kglob(3:6,3:6)+Kloc2;


% % Load vector (in a global coordinate systems)
fglob=zeros(6,1);
fglob(5)=F;

% Boundary conditions (in a global coordinate systems)
% v1=0,theta1=0 
% Remaining (free) dofs (nodal displacements) are related to indeces 3,4,5,6 (1,2 (v1,theta1) are fixed)
% Remanininf DOFs of system = 6 -2 (number of constraints) = 4
Kglobred=zeros(4,4);
fglobred=zeros(4,1);
Kglobred=Kglob([3:6],[3:6]);
fglobred=fglob([3:6]);

% Let's solve determinant
det(Kglobred)

% Let's solve displacements uglob 
uglobred=Kglobred\fglobred

% Pre processing
% Let's gather all displacements (solved and fixed)
ugloball=zeros(4,1);
ugloball([3:6])=uglobred

% Let's gather elemental displacement vectors
ug1=ugloball(1:4);
ug2=ugloball(3:6);

uloc1=ug1;
uloc2=ug2;

% Let's solve member forces (in a local coordinate system)
Floc1=Kloc1*uloc1
Floc2=Kloc2*uloc2

% Analytical solution based on the Euler-Bernoulli beam theory (deflection)
vmaxAnalytic=F*L^3/(3*E*Iz)

% Let's draw deflection curve. Now x is at free end.
x=L:-L/100:0;

for ii=1:length(x)
    vAnalytic(ii)=F*L^3/(6*E*Iz)*(2-3*x(ii)/L+x(ii)^3/L^3);
end

figure(1)
plot(x,vAnalytic,'b-')
hold on
plot([0,L/2,L], ugloball([5,3,1]),'r*')


