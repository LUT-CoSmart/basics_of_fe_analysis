% Basics of Finite Element Method
% Marko Matikainen
% Example of differential equation for a tapered rod 
% and MATLAB symbolic toolbox

% Introduce variablses
syms u(x) A(x) E F L A1 A2 r1 r2

% Linearily (cross section variates) tapered rod. You can find equation for A(x) by yourself.
%A(x) = A1*(1 - x/L) + A2*x/L;       % but here it is given.

% Radius varies linearly along the rod
r = r1*(1 - x/L) + r2*x/L;

% Circular cross-sectional area
A = pi*r^2;

epsilon = diff(u);                  % strain epsilon = du/dx                    
%DE = -diff(E*A(x)*epsilon) == 0     % Differential equation
%bc = [u(0) == 0; E*A(L)*epsilon(L) == F] % Boundary conditions

% Differential equation
DE = -diff(E*A*epsilon,x) == 0;

% Boundary conditions
% Fixed at x = 0 and axial force F at x = L
bc = [u(0) == 0;
      subs(E*A*epsilon,x,L) == F];


u = dsolve(DE,bc)                   % Solves DE with given bc
u = expand(u)                       
u = simplify(u)
latex(u)