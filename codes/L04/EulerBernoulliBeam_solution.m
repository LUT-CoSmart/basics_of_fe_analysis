% Code derives solution for E-B eq.

syms x L E I q0
syms v(x)

% Euler-Bernoulli beam equation
ode = E*I*diff(v,x,4) == q0;

% Boundary conditions for a cantilever beam
bc1 = v(0) == 0;                         % displacement
bc2 = subs(diff(v,x),x,0) == 0;         % rotation
bc3 = subs(diff(v,x,2),x,L) == 0;       % bending moment
bc4 = subs(diff(v,x,3),x,L) == 0;       % shear force

% Analytical solution
vanal = simplify(dsolve(ode,[bc1 bc2 bc3 bc4]))