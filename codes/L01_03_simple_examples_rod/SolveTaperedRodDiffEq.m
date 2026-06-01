% Basics of Finite Element Method
% Marko Matikainen
% Example of differential equation for a tapered rod 
% and MATLAB symbolic toolbox

% Introduce variablses
syms u(x) A(x) E F L A1 A2

% Linearily tapered rod. You can find equation for A(x) by yourself.
A(x) = A1*(1 - x/L) + A2*x/L;       % but here it is given.

epsilon = diff(u);                  % strain epsilon = du/dx                    
DE = -diff(E*A(x)*epsilon) == 0     % Differential equation
bc = [u(0) == 0; E*A(L)*epsilon(L) == F] % Boundary conditions

u = dsolve(DE,bc)                   % Solves DE with given bc
u = expand(u)                       
u = simplify(u)
latex(u)