% BK10A6400 Basics of FE Analysis (FEMBasics2025autumn)
% Teacher in charge: Marko Matikainen (MKM)
% Example: Solving the displacement of a tapered rod structure with four rod elements
%
% Coded by MKM for student's usage in the FEMBasics2025 course

%             x
%           |--> 
%
%          /u1      u2      u3      u4      u5
%          /o-------o-------o-------o-------o-> F
%          /N1  E1  N2  E2  N3  E3  N4  E4  N5
%
%           |<------------- L ------------->|  
%
%   Boundary conditions:
%       u1(x=0) = 0
%       f(x=L) = F
%
% Modified by ChatGPT 2.6.2026 (not verified!!)


clear all;
clc;
close all;
format shortG;

% Initial parameters
E = 210000;      % Young's modulus [MPa] = [N/mm^2] 
r1 = 50;         % Radius at x=0, fixed end [mm]
r2 = 25;         % Radius at x=L, free end [mm]
   
A1 = pi*r1^2;    % Cross-sectional area at x=0, fixed end [mm^2]
A2 = pi*r2^2;    % Cross-sectional area at x=L, free end [mm^2]

L = 1000;        % Length of the structure [mm]
F = 1000;        % Applied force [N]

% Number of elements and nodes
nel = 4;
nnod = nel + 1;

% Element length
Le = L/nel;

% Calculate the average cross-sectional area for each element
% The cross-sectional areas are evaluated at the element midpoints
A1avg = LinCrossSec(A1,A2,L,1/8*L);
A2avg = LinCrossSec(A1,A2,L,3/8*L);
A3avg = LinCrossSec(A1,A2,L,5/8*L);
A4avg = LinCrossSec(A1,A2,L,7/8*L);

% Element stiffness matrices
Kloc1 = KlocRod(E,A1avg,Le);
Kloc2 = KlocRod(E,A2avg,Le);
Kloc3 = KlocRod(E,A3avg,Le);
Kloc4 = KlocRod(E,A4avg,Le);

% Transformation, if needed
% The elements are parallel to the global coordinate system, so no transformation is needed
Kglob1 = Kloc1;
Kglob2 = Kloc2;
Kglob3 = Kloc3;
Kglob4 = Kloc4;

% Initialise the global stiffness matrix and force vector
Kglob = zeros(nnod,nnod);   % 5 nodes, 1 DOF per node -> 5 DOFs
fglob = zeros(nnod,1);

% Apply force at the free end
fglob(5) = F;

% Assemble the global stiffness matrix
Kglob(1:2,1:2) = Kglob(1:2,1:2) + Kglob1;
Kglob(2:3,2:3) = Kglob(2:3,2:3) + Kglob2;
Kglob(3:4,3:4) = Kglob(3:4,3:4) + Kglob3;
Kglob(4:5,4:5) = Kglob(4:5,4:5) + Kglob4;

% Apply boundary condition at the fixed end
% u1 = 0, so the first row and first column are removed
Kglobc = Kglob(2:5,2:5);
fglobc = fglob(2:5);

% Solve the nodal displacements
uglobc = Kglobc\fglobc;

% Full displacement vector including the fixed displacement u1 = 0
uglob = zeros(nnod,1);
uglob(2:5) = uglobc;

disp('Global stiffness matrix Kglob:')
Kglob

disp('Reduced stiffness matrix Kglobc:')
Kglobc

disp('Global nodal displacement vector uglob:')
uglob

% Analytical solution at the free end based on linear theory
if A1 == A2    
    umax = F*L/(E*A1);  
else
    umax = (F*L/(E*(A2 - A1))) * log(A2/A1);
end

disp('Analytical displacement at the free end:')
umax

disp('FE displacement at the free end:')
u_FE_free_end = uglob(end)

disp('Absolute error at the free end:')
error_abs = abs(umax - u_FE_free_end)

disp('Relative error [%]:')
error_rel = abs(umax - u_FE_free_end)/abs(umax)*100