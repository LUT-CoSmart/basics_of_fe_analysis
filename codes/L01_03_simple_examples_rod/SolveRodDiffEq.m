% % BK10A6400 Basics of FE Analysis (FEMBasics2025)
% Teacher in charge: Marko Matikainen 
% The code solves a rod DE 
% Coded by MKM for student's usage in the FEMBasics2025 course

% Theory recap: DE for a rod
% https://appliedmechanics.ju.se/TheRodEquation/


syms rho g u(x) E A L F

q = rho*g;
du = diff(u);                       % strain epsilon = du/dx
DE = -diff(E*A*du) == q             % Differential equation
bc = [u(0) == 0; E*A*du(L) == F]    % Boundary conditions

u = dsolve(DE,bc)                   % Solves DE with given bc
u = expand(u)

