% BK10A6400 Basics of FE Analysis (FEMBasics2022)
% Teacher in charge: Marko Matikainen (MKM)
% Solution for the Bonus Task 3.
%
% Goal: The code solves global displacements at node 2 and member forces of a simple beam structure. 
% Coded by MKM for student's usage in the FEMBasics2022 course

%          N1     N2       N3
%          L1,I1     L2,I2
%          \      |F
%          \o-----o--------o
%          \  (1)     (2) o^o  
%
%  v1,theta1   v2,theta2  v3,theta3

% So-called connectivity matrix
% Ele1 N1 N2   v1,theta1,v2,theta2
% Ele2 N2 N3   v2,theta2,v3,theta3
% so node 2 is common node (and then displacements v2,theta2) for the both beam elements

% First matlab-based FE code for a simple beam structure
clear all;
close all;
clc;
%format long;
format shortG;

% Units are mm and N and radians
% Beam 1
E1=210000;
I1=4762;
L1=1000;
alpha1=0;

% Beam 2
E2=210000;
I2=9524;
L2=2000;
alpha2=0;

% Force (in vertical (global) direction at node 2)
F=-100;

% Elemental stiffness matrices
Kloc1=KlocEB4DOF(E1,I1,L1);
Kloc2=KlocEB4DOF(E2,I2,L2);

% Transformation into a global coordinate system
% Transformation matrix for a beam 1 
T1=TEB4DOFs(alpha1);
% Transformation matrix for a beam 2
T2=TEB4DOFs(alpha2);

% Global elemental stiffness matrix for a rod 1
Kglob1=T1'*Kloc1*T1;
% Global elemental stiffness matrix for a rod 2
Kglob2=T2'*Kloc2*T2;

% As you can see, local stiffness matrices have zero rows and columns due
% to transformation to global coordinate system and paraller coordinate
% systems so no point to make transformation because zero rows and columns 
% need to eliminated. Note that this is just special case. For 6 DOFs beam
% element   later, use transformation.

%Let's use Klocs directly:
Kglob=zeros(6,6);
Kglob(1:4,1:4)=Kloc1;
Kglob(3:6,3:6)=Kglob(3:6,3:6)+Kloc2;

% % If this previous assembling way in matlab is somehow difficult,
% % you can always do:
% Kglob(1:4,1:4)=Kloc1;
% Kglob(3:6,3:6)=Kloc2;
% % but now stiffnesses of a common node is not computed properly
% % So lets do it explicitely
% Kglob(3,3)=Kloc1(3,3)+Kloc2(1,1);
% Kglob(3,4)=Kloc1(3,4)+Kloc2(1,2);
% Kglob(4,3)=Kloc1(4,3)+Kloc2(2,1);
% Kglob(4,4)=Kloc1(4,4)+Kloc2(2,2);

% % Load vector (in a global coordinate systems)
fglob=zeros(6,1);
fglob(3)=F;

% Boundary conditions (in a global coordinate systems)
% v1=0,theta1=0,v3=0 
% Remaining (free) dofs (nodal displacements) are related to indeces 3,4,6 (1,2,5 (v1,theta1,v3) are fixed)
% Remanininf DOFs of system = 6 -3 (number of constraints) = 3
Kglobred=zeros(3,3);
fglobred=zeros(3,1);
Kglobred=Kglob([3:4,6],[3:4,6]);
fglobred=fglob([3:4,6]);

% Let's solve determinant
det(Kglobred)

% Let's solve displacements uglob 
uglobred=Kglobred\fglobred

% Next lines are just for see why it is beneficial to use
% backslash for solvig a system of linear equations uglobred=Kglobred\fglobred 
% Comment lines out and see computations time for these different opearations

% tic  % tic - starts clock, toc - ends clock 
% uglobred=inv(Kglobred)*fglobred
% toc
% 
% tic
% uglobred=Kglobred^(-1)*fglobred
% toc
% 
% tic
% uglobred=Kglobred\fglobred
% toc

% Let's gather all displacements (solved and fixed)
ugloball=zeros(6,1);
ugloball([3:4,6])=uglobred

% Let's gather elemental displacement vectors
ug1=ugloball(1:4);
ug2=ugloball(3:6);

% % Now no need to make transformation
% uloc1=T1*ug1;
% uloc2=T2*ug2;

uloc1=ug1
uloc2=ug2

% Let's solve member forces (in a local coordinate system)
Floc1=Kloc1*uloc1
Floc2=Kloc2*uloc2


% Post processing
% ###############################################################
% ###############################################################
%% Simple visualization using only nodal displacements

scaleFactor = 20;   % visualization scale factor

% Nodal x-coordinates
xnod = [0, L1, L1+L2];

% Nodal vertical displacements
vnod = [ugloball(1), ugloball(3), ugloball(5)];

figure(1)
hold on

% Undeformed beam
plot(xnod,[0 0 0],'k--','LineWidth',1.5)

% Deformed beam using only nodal displacements
plot(xnod,scaleFactor*vnod,'-o','LineWidth',2)


grid on
box on
axis equal

xlabel('x [mm]')
ylabel('Scaled vertical displacement')

title('Finite element solution')

legend('Undeformed beam', ...
       'Deformed beam', ...
       'Location','best')

fprintf('Visualization scale factor = %g\n',scaleFactor)



%% Visualization using shape functions

scaleFactor = 50;
nplot = 50;
axis equal

% ----------------------------------
% Element 1
% ----------------------------------

x1loc = linspace(0,L1,nplot);
vplot1 = zeros(size(x1loc));

for ii = 1:length(x1loc)
    N = Shapef_EB4DOF(x1loc(ii),L1);
    vplot1(ii) = N*uloc1;
end

x1glob = x1loc;


% ----------------------------------
% Element 2
% ----------------------------------

x2loc = linspace(0,L2,nplot);
vplot2 = zeros(size(x2loc));

for ii = 1:length(x2loc)
    N = Shapef_EB4DOF(x2loc(ii),L2);
    vplot2(ii) = N*uloc2;
end

x2glob = L1 + x2loc;


% ----------------------------------
% Plot
% ----------------------------------

figure(2)
hold on

% Undeformed beam
plot([0 L1 L1+L2],[0 0 0],'k--')

% FE displacement field
plot(x1glob,scaleFactor*vplot1,'LineWidth',2)
plot(x2glob,scaleFactor*vplot2,'LineWidth',2)

grid on
axis equal

xlabel('x [mm]')
ylabel('Scaled displacement [mm]')

legend('Undeformed',...
       'Element 1',...
       'Element 2')