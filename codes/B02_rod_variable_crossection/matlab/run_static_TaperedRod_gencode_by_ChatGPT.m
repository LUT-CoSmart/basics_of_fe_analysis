% BK10A6400 Basics of FE Analysis (FEMBasics2025autumn)
% Teacher in charge: Marko Matikainen (MKM)
% Example: Solving the displacement of a tapered rod structure
% with an arbitrary number of rod elements
%
% Coded by MKM for student's usage in the FEMBasics2025 course

clear all;
clc;
close all;
format shortG;

% Initial parameters
E  = 210000;      % Young's modulus [MPa] = [N/mm^2]
r1 = 50;          % Radius at x=0, fixed end [mm]
r2 = 25;          % Radius at x=L, free end [mm]

A1 = pi*r1^2;     % Cross-sectional area at x=0 [mm^2]
A2 = pi*r2^2;     % Cross-sectional area at x=L [mm^2]

L = 1000;         % Length of the structure [mm]
F = 1000;         % Applied force [N]

% Number of elements
nel = 128;          % Change this value: 2, 4, 8, 16, ...

% Number of nodes
nnod = nel + 1;

% Element length
Le = L/nel;

% Initialisation of global stiffness matrix and force vector
Kglob = zeros(nnod,nnod);   % 1 DOF per node
fglob = zeros(nnod,1);

% Applied force at the free end
fglob(end) = F;

% Assembly of the global stiffness matrix
for e = 1:nel

    % Element midpoint coordinate
    xmid = (e - 0.5)*Le;

    % Average cross-sectional area of the element
    Aavg = LinCrossSec(A1,A2,L,xmid);

    % Element stiffness matrix
    Kloc = KlocRod(E,Aavg,Le);

    % Element degrees of freedom
    edof = [e e+1];

    % Assembly into the global stiffness matrix
    Kglob(edof,edof) = Kglob(edof,edof) + Kloc;

end

% Boundary condition at the fixed end
% u1 = 0, therefore the first row and first column are removed
freeDofs = 2:nnod;

Kglobc = Kglob(freeDofs,freeDofs);
fglobc = fglob(freeDofs);

% Solve unknown displacements
uglobc = Kglobc\fglobc;

% Full displacement vector including fixed displacement u1 = 0
uglob = zeros(nnod,1);
uglob(freeDofs) = uglobc;

% Analytical solution at the free end
if A1 == A2
    umax = F*L/(E*A1);
else
    umax = (F*L/(E*(A2 - A1))) * log(A2/A1);
end

% FE displacement at the free end
u_FE_free_end = uglob(end);

% Error at the free end
error_abs = abs(umax - u_FE_free_end);
error_rel = error_abs/abs(umax)*100;

% Print results
disp('Global stiffness matrix Kglob:')
Kglob

disp('Global force vector fglob:')
fglob

disp('Global displacement vector uglob:')
uglob

disp('Analytical displacement at the free end:')
umax

disp('FE displacement at the free end:')
u_FE_free_end

disp('Absolute error at the free end:')
error_abs

disp('Relative error [%]:')
error_rel


% -------------------------------------------------------------------------
% Visualisation of FE and analytical solution
% -------------------------------------------------------------------------

% Nodal coordinates for FE solution
xnod = linspace(0,L,nnod);

% Dense x-vector for analytical solution
x_analytical = linspace(0,L,200);

% Analytical displacement as a function of x
if A1 == A2
    u_analytical = F*x_analytical/(E*A1);
else
    A_x = A1 + (A2 - A1)*x_analytical/L;
    u_analytical = (F*L/(E*(A2 - A1))) * log(A_x/A1);
end

figure;
plot(xnod,uglob,'o-','LineWidth',1.5);
hold on;
plot(x_analytical,u_analytical,'--','LineWidth',1.5);
grid on;

xlabel('x [mm]');
ylabel('Displacement u [mm]');
title(['Axial displacement, number of elements = ',num2str(nel)]);

legend('FE solution','Analytical solution','Location','northwest');


% -------------------------------------------------------------------------
% Local functions
% These can also be saved as separate .m files if preferred
% -------------------------------------------------------------------------

function A = LinCrossSec(A1,A2,L,x)
% Linear interpolation of cross-sectional area

A = A1 + (A2 - A1)*x/L;

end


function Kloc = KlocRod(E,A,L)
% Local stiffness matrix of a two-node rod element

Kloc = E*A/L*[ 1 -1;
              -1  1 ];

end