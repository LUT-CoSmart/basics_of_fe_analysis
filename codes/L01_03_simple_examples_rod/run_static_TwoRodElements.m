% BK10A6400 Basics of FE Analysis (FEMBasics2025autumn)
% Teacher in charge: Marko Matikainen (MKM)
% Example: Solving displacement of a rod structure with two rod elements
%
% Goal: The code solves global displacements at node 2 and member forces of a simple rod structure. 
% Coded by MKM for student's usage in the FEMBasics2025 course

%             x
%           |--> 
%
%            L1,A1,E1   L2,A2,E2 
%          /u1       u2       u3 
%          /o--------o--------o-> F
%          /N1  E1   N2  E2   N3  
%
%
%           |<-------L------->|  
%
%   bc: u1(x=0)=0, 
%       f(x=L)=F
%

clear all;
clc;
close all;
format shortG;

% Grazy initials
% E=1;A=1;L=1;F=1;


% Something that makes sense
E=210000;       % Young's modulus [Mpa] 
D=50;           % radius [mm]
A=pi*D^2/4;     % Cross section area [mm^2]
L=1000;         % Length of structure [mm]
F=1000;         % Applied force [N]

% Element 1
L1=L/2;
A1=A;
E1=E;
% Element 2
L2=L/2;
A2=A;
E2=E;

Kloc1 = KlocRod(E1,A1,L1);
Kloc2 = KlocRod(E2,A2,L2);

% Transformation if needed
% Elements are paraller with the global coordinate system...no need.
Kglob1=Kloc1;
Kglob2=Kloc2;

% It is possible to make transformation but then, zero rows and columns
% need to be eliminated away. Check how it goes...

% Global FE vector/matrices for assembling
% Intialistion
Kglob=zeros(3,3);       % 3 nodes, 1 DOF per node -> 3 DOFS
fglob=zeros(3,1);

fglob(3)=F;

% Assembling
% Option 1
%Kglob(1:2,1:2)=Kloc1;
%Kglob(2:3,2:3)=Kloc2;
%Kglob(2,2)=Kloc1(2,2)+Kloc2(1,1);

% Option 2
Kglob(1:2,1:2)=Kglob1;      % Substitute element 1 in a global system stiffness matrix
Kglob(2:3,2:3)=Kglob(2:3,2:3)+Kglob2; % substitute element 2 and sum stiffnesses from elems 1 and 2

% Boundary conditions for clambed end
Kglobc=Kglob(2:3,2:3);
fglobc=fglob(2:3);

% Diplacement
uglobc=Kglobc\fglobc        % The most efficient way to solve the system of linear equations in MATLAB. It is recommended by MATLAB

% uglobc=Kglobc^-1*fglobc   % this works also but use backslash instead
% uglobc=inv(Kglobc)*fglobc % this works also but use backslash instead

% Analytical solution based on linear theory at free end.
% Jus derive it from sigma=F/A, epsilon=deltau/L
deltau=F*L/(E*A);




