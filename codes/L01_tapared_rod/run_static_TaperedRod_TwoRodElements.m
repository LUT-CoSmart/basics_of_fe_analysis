% BK10A6400 Basics of FE Analysis (FEMBasics2024autumn)
% Teacher in charge: Marko Matikainen (MKM)
% Example: Solving displacement of a tapered rod structure with two rod elements
%
% Coded by MKM for student's usage in the FEMBasics2024 course

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

% Initial parameters
E=210000;       % Young's modulus [Mpa] 
r1=50;          % radius at x=0 (fixed end) [mm]
r2=25;          % radius at x=L (free end) [mm]
   
A1=pi*r1^2;     % Cross section area at x=0 (fixed end)
A2=pi*r2^2;     % Cross section area at x=L (free end)

L=1000;         % Length of structure [mm]
F=1000;         % Applied force [N]

% Calculates elemental average cross section
A1avg = LinCrossSec(A1,A2,L,L/4);       % Averaged cross section at x=L/4 
A2avg = LinCrossSec(A1,A2,L,3/4*L);     % Averaged cross section at x=3/4*L

% Elemental stiffness matrices
Kloc1 = KlocRod(E,A1avg,L/2);
Kloc2 = KlocRod(E,A2avg,L/2);

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

% Applied force at free end
fglob(3)=F;

% Assembling
% Option 1
%Kglob(1:2,1:2)=Kloc1;
%Kglob(2:3,2:3)=Kloc2;
%Kglob(2,2)=Kloc1(2,2)+Kloc2(1,1);

% Option 2
Kglob(1:2,1:2)=Kglob1;      % Substitute element 1 in a global system stiffness matrix
Kglob(2:3,2:3)=Kglob(2:3,2:3)+Kglob2; % substitute element 2 and sum stiffnesses from elems 1 and 2

% Boundary conditions for clambed (fixed) end
Kglobc=Kglob(2:3,2:3);
fglobc=fglob(2:3);

% Diplacement
uglobc=Kglobc\fglobc        % The most efficient way to solve the system of linear equations in MATLAB. It is recommended by MATLAB

% uglobc=Kglobc^-1*fglobc   % this works also but use backslash instead
% uglobc=inv(Kglobc)*fglobc % this works also but use backslash instead


% Analytical solution based on linear theory at free end.
if A1==A2,    
    umax=F*L/(E*A1);  
else
    umax = (F * L / (E * (A2 - A1))) * log(A2 / A1);
end

% Analytical solution as function of x
%x=L;
%umax = (F*L*(log(-A1*L) - log(A1*x - A2*x - A1*L)))/(E*(A1 - A2));

umax