% BK10A6400 Basics of FE Analysis (FEMBasics2025)
% Teacher in charge: Marko Matikainen (MKM)
% Solution for the Bonus Task 4.
%
% Goal: The code solves global displacements at node 2 and member forces of a simple beam structure. 
% Coded by MKM for student's usage in the FEMBasics2022 course

%          N1     N2       N3
%          L1,I1     L2,I2
%          \      |F
%          \o-----o--------o
%          \  (1)     (2) o^o  
%
%  v1,theta1   v2,theta2  v3,theta3

% So-called connectivity matrix
% Ele1 N1 N2   v1,theta1,v2,theta2
% Ele2 N2 N3   v2,theta2,v3,theta3
% so node 2 is common node (and then displacements v2,theta2) for the both beam elements

% First matlab-based FE code for a simple beam structure
clear all;
close all;
clc;
format long;

% Units are mm and N and radians
% Beam 1
E1=210000;
I1=4762;
L1=1000;
alpha1=0;

% Beam 2
E2=210000;
I2=9524;
L2=2000;
alpha2=0;

% Force (in vertical (global) direction at node 2)
F=-100;

% Elemental stiffness matrices
Kloc1=KlocEB4DOF(E1,I1,L1);
Kloc2=KlocEB4DOF(E2,I2,L2);

% No need for transformatio so let's use Klocs directly. The system has
% three nodes and then the system's global stiffness matrix is a size of
% 6x6
Kglob=zeros(6,6);
%Kglob(?,?)=?       % TIP: make assembling. Take a look at the lecture slodes and previous examples.
%Kglob(?,?)=?       


% % Load vector (in a global coordinate systems)
fglob=zeros(6,1);       % TIP: add plied force in a correct location of a vector of external forces.
% fglob(?)=?;

% Boundary conditions (in a global coordinate systems)
% v1=0,theta1=0,v3=0 
% Remaining (free) dofs (nodal displacements) are related to indeces 3,4,6 (1,2,5 (v1,theta1,v3) are fixed)
% Remanininf DOFs of system = 6 -3 (number of constraints) = 3
Kglobred=zeros(3,3);
fglobred=zeros(3,1);
% Kglobred=Kglob([?],[?]);        % TIP: take boundary constraints into account.
% fglobred=fglob([?]);

% Let's solve determinant.
det(Kglobred)

% Let's solve displacements uglob 
uglobred=Kglobred\fglobred

% Let's gather all displacements (solved and fixed)
ugloball=zeros(6,1);
% ugloball([?])=uglobred;       % TIP: add solved displacements uglobred
% into correct locations

% Let's gather elemental displacement vectors
% ug1=ugloball(?);              % TIP: pick up correct stiffnesses
% ug2=ugloball(?);

uloc1=ug1;
uloc2=ug2;

% Let's solve member forces (in a local coordinate system)
% Floc1=                        % TIP : Solve member forces, take a look at the
% lecture notes
% Floc2=


