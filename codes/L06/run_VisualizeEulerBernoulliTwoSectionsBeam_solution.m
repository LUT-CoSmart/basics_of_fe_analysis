%% Analytical solution for Bonus Task 3
%
% Uses analytical Euler-Bernoulli solutions:
%
%   v1anal.m      -> deflection of beam section 1
%   v2anal.m      -> deflection of beam section 2
%   theta1anal.m  -> rotation of beam section 1
%   theta2anal.m  -> rotation of beam section 2
%
% The external force F acts at x = L1.

clear
clc
close all


%% ------------------------------------------------------------------------
% Parameters of the Bonus Task
% -------------------------------------------------------------------------

E  = 210000;       % Young's modulus [MPa = N/mm^2]

I1 = 4762;         % Second moment of area, beam section 1 [mm^4]
I2 = 9524;         % Second moment of area, beam section 2 [mm^4]

L1 = 1000;         % Length of beam section 1 [mm]
L2 = 2000;         % Length of beam section 2 [mm]

F  = 100;          % External force [N]


%% ------------------------------------------------------------------------
% Vertical displacement at the location of force F
%
% Force F acts at:
%
%       x = L1
%
% The displacement can be evaluated using either v1anal or v2anal
% because displacement is continuous at x = L1.
% -------------------------------------------------------------------------

uY_v1 = v1anal(L1,E,I1,I2,L1,L2,F);
uY_v2 = v2anal(L1,E,I1,I2,L1,L2,F);


%% ------------------------------------------------------------------------
% Rotation at the location of force F
%
% Rotation:
%
%       theta_Z = dv/dx
%
% The rotation can be evaluated using either theta1anal or theta2anal
% because rotation is continuous at x = L1.
% -------------------------------------------------------------------------

thetaZ_v1 = theta1anal(L1,E,I1,I2,L1,L2,F);
thetaZ_v2 = theta2anal(L1,E,I1,I2,L1,L2,F);


%% ------------------------------------------------------------------------
% Print analytical solution at the loaded node
% -------------------------------------------------------------------------

fprintf('\n')
fprintf('ANALYTICAL EULER-BERNOULLI SOLUTION\n')
fprintf('------------------------------------\n')

fprintf('Displacement from beam section 1: uY     = %12.8f mm\n',uY_v1)
fprintf('Displacement from beam section 2: uY     = %12.8f mm\n',uY_v2)

fprintf('Rotation from beam section 1:     thetaZ = %12.8f rad\n',thetaZ_v1)
fprintf('Rotation from beam section 2:     thetaZ = %12.8f rad\n',thetaZ_v2)


%% ------------------------------------------------------------------------
% Check continuity at x = L1
% -------------------------------------------------------------------------

fprintf('\n')
fprintf('Continuity checks at x = L1\n')
fprintf('---------------------------\n')

fprintf('Difference in displacement = %e mm\n', ...
        abs(uY_v1-uY_v2))

fprintf('Difference in rotation     = %e rad\n', ...
        abs(thetaZ_v1-thetaZ_v2))


%% ------------------------------------------------------------------------
% Coordinates for plotting the analytical deflection curve
% -------------------------------------------------------------------------

nplot = 200;

%x1 = linspace(0,L1,nplot);
%x2 = linspace(L1,L1+L2,nplot);

% lets presend defelctiosn outside of physical area:
x1 = linspace(-L1,L1+L2,nplot);
x2 = linspace(0,L1+L2,nplot);


%% ------------------------------------------------------------------------
% Evaluate analytical deflection curve
% -------------------------------------------------------------------------

v1 = v1anal(x1,E,I1,I2,L1,L2,F);
v2 = v2anal(x2,E,I1,I2,L1,L2,F);

%% Stationary points

x_v1 = x_v1max(I1,I2,L1,L2);
x_v2 = x_v2max(I1,I2,L1,L2);

disp('Stationary points for beam section 1:')
disp(x_v1)

disp('Stationary points for beam section 2:')
disp(x_v2)


%% Check which stationary points belong to the beam sections

% Beam section 1:
% 0 <= x <= L1

% Beam section 2:
% L1 <= x <= L1 + L2


%% For this problem the valid stationary point is in beam section 2

x_max = x_v2(2);
uY_max = v2anal(x_max,E,I1,I2,L1,L2,F);


fprintf('\n')
fprintf('Maximum vertical displacement:\n')
fprintf('x      = %12.4f mm\n',x_max)
fprintf('uY_max = %12.8f mm\n',uY_max)



%% ------------------------------------------------------------------------
% Plot analytical deflection curve
% -------------------------------------------------------------------------

figure

plot(x1,v1,'LineWidth',1.5)
hold on

plot(x2,v2,'LineWidth',1.5)

plot(L1,uY_v1,'o','MarkerSize',7)

grid on
box on

xlabel('x [mm]')
ylabel('u_Y(x) [mm]')

title('Analytical Euler-Bernoulli beam solution')

legend('Beam section 1', ...
       'Beam section 2', ...
       'Location of F', ...
       'Location','best')


%% ------------------------------------------------------------------------
% Final result at the loaded node
% -------------------------------------------------------------------------

fprintf('\n')
fprintf('RESULT AT THE LOADED NODE\n')
fprintf('-------------------------\n')

fprintf('uY     = %12.8f mm\n',uY_v1)
fprintf('thetaZ = %12.8f rad\n',thetaZ_v1)

fprintf('\n')