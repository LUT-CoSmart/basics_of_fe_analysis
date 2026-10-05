% Analytical solution for two-section E-B beam
% Point force F acts at x = L1
% Check the solution!

clear
clc

syms x L1 L2 E I1 I2 F
syms v1(x) v2(x)

% Beam equations: no distributed load q(x)
ode1 = E*I1*diff(v1,x,4) == 0;
ode2 = E*I2*diff(v2,x,4) == 0;


%% Boundary conditions

% Left end: clamped
bc1 = v1(0) == 0;
bc2 = subs(diff(v1,x),x,0) == 0;

% Right end: hinged
bc3 = v2(L1+L2) == 0;
bc4 = subs(diff(v2,x,2),x,L1+L2) == 0;


%% Conditions at x = L1

% Displacement continuity
bc5 = v1(L1) == v2(L1);

% Rotation continuity
bc6 = subs(diff(v1,x),x,L1) == ...
      subs(diff(v2,x),x,L1);

% Bending moment continuity
bc7 = E*I1*subs(diff(v1,x,2),x,L1) == ...
      E*I2*subs(diff(v2,x,2),x,L1);

% Jump in shear force due to point force F
bc8 = E*I2*subs(diff(v2,x,3),x,L1) - ...
      E*I1*subs(diff(v1,x,3),x,L1) == -F;


%% Analytical solution

[v1anal,v2anal] = dsolve( ...
    [ode1 ode2], ...
    [bc1 bc2 bc3 bc4 bc5 bc6 bc7 bc8]);

v1anal = simplify(v1anal);
v2anal = simplify(v2anal);


% Save analytical deflection solutions as MATLAB functions
matlabFunction(v1anal,'file','v1anal','vars',{x,E,I1,I2,L1,L2,F});
matlabFunction(v2anal,'file','v2anal','vars',{x,E,I1,I2,L1,L2,F});

% Rotations from deflections (small angle assumption)
theta1anal = simplify(diff(v1anal,x));
theta2anal = simplify(diff(v2anal,x));

matlabFunction(theta1anal,'file','theta1anal','vars',{x,E,I1,I2,L1,L2,F});
matlabFunction(theta2anal,'file','theta2anal','vars',{x,E,I1,I2,L1,L2,F});

% Let's find maximum deflection for each section
% d v1 / dx = 0 -> x = 
dv1dx = simplify(diff(v1anal,x));   % Note, this is same than rotation /angle 
dv2dx = simplify(diff(v2anal,x));    % Note, this is same than rotation /angle 
x_v1max = solve(dv1dx == 0,x);      % statinary points
x_v2max = solve(dv2dx == 0,x);     % statinary points 

% Save stationary-point solutions as MATLAB functions
matlabFunction(x_v1max,'file','x_v1max','vars',{I1,I2,L1,L2});
matlabFunction(x_v2max,'file','x_v2max','vars',{I1,I2,L1,L2});


