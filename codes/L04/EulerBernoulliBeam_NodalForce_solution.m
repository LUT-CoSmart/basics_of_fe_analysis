% Code derives solution for E-B eq.
% for cantilever beam under applied force F

syms x L E I F
syms v(x)

% Beam equation: no distributed load q(x)
ode = E*I*diff(v,x,4) == 0;

% Boundary conditions
bc1 = v(0) == 0;
bc2 = subs(diff(v,x),x,0) == 0;
bc3 = subs(diff(v,x,2),x,L) == 0;
bc4 = E*I*subs(diff(v,x,3),x,L) == F;

vanal = simplify(dsolve(ode,[bc1 bc2 bc3 bc4]))