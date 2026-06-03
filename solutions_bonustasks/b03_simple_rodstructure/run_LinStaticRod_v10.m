% BK10A6400 Basics of FE Analysis (FEMBasics2023autumn)
% Teacher in charge: Marko Matikainen (MKM)
% Final solution code for the Bonus Task 2.
%
% Goal: The code solves global displacements at node 2 and member forces of a simple rod structure. 
% Coded by MKM for student's usage in the FEMBasics2023 course

%           o N1
%            \
%       Elem1 \   
%              o  N2  
%       Elem2 /|
%            /  F 
%           o N3


% First matlab-based FE code for a simple truss structure
clear all;
close all;
clc;
format shortG;  % just to display values in more convient manner

% Units are mm and N and radians
% Rod 1
E1=210000;
A1=100;
L1=1000;
alpha1=-45*pi/180;
%alpha1=315*pi/180;

% Rod 2
E2=210000;
A2=50;
L2=1000;
alpha2=-135*pi/180;

% Force (in vertical (global) direction at node 2)
F=-100000;


% Transformation into a global coordinate system
% Transformation matrix for a rod 1 
T1=TRod(alpha1);
% Transformation matrix for a rod 2
T2=TRod(alpha2);

% Local elemental stiffness matrices
Kloc1=KlocRod(E1,A1,L1);
Kloc2=KlocRod(E2,A2,L2);

% Global elemental stiffness matrix for a rod 1
Kglob1=T1'*Kloc1*T1;
% Global elemental stiffness matrix for a rod 2
Kglob2=T2'*Kloc2*T2;

% Global stiffness matrix assembling for a whole structure (6x6)
Kglob=zeros(6,6);           % Initialize zero matrix
Kglob(1:4,1:4)=Kglob1;      % Substitute element 1 in a global system stiffness matrix
Kglob(3:6,3:6)=Kglob(3:6,3:6)+Kglob2; % substitute element 2 and sum stiffnesses from elems 1 and 2

% Load vector (in a global coordinate systems)
fglob=zeros(6,1);           % Initialization
fglob(4)=F;                 % Substitute nodal force with respectt to a correct DOF

% Boundary conditions (in a global coordinate systems)
% After bc's, global stiffness matrix is 2x2 and external force vector is
% 2x1.
% Boundary conditions u1,v1,u3,v3 are fixed so corresponding rows and
% columns are eliminated from Kglob and fglob. 
Kglobc=Kglob(3:4,3:4);      % just substitute correct rows and columns
fglobc=fglob(3:4);

% Now it is good to check that the stiffness matrix has full rank and/or det
% ~= 0 to be sure that it is invertible.
rank(Kglobc);
det(Kglobc);

% Solve nodal displacement (in a global coordinate systems)
uglobc=Kglobc\fglobc;

% Let's collect all nodal displacements (solved and fixed)
ugloball=zeros(6,1);
ugloball(3:4)=uglobc;

% Makes elemental nodal displacement vectors
uglob1=ugloball(1:4);
uglob2=ugloball(3:6);

% Solve element's axial displacements (in a local coordinate systems)
uloc1=T1*uglob1;
uloc2=T2*uglob2;

% Let's solve member forces f=K u
floc1=Kloc1*uloc1;
floc2=Kloc2*uloc2;

% Next ones are asked to return in Moodle
disp('Determinant of Kglob')
det(Kglob)             % Determinant of unconstrained global stiff matrix (should be 0 or really close)
disp('Determinant of Kglobc')
det(Kglobc)            % Determinant of constrained global stiff matrix (>>0)
disp('Displacement at node 2')
ugloball(3:4)          % Displacements at node 2
disp('Local elemental displacements')
uloc1
uloc2
disp('Member forces')
floc1
floc2






