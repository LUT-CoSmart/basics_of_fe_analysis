% % BK10A6400 Basics of FE Analysis (FEMBasics2025autumn)
% Teacher in charge: Marko Matikainen (MKM)
% Template for the Bonus Task 4.

% Let's test developed functions for the internal forces.

clear all;
close all;
clc;

% Initial values
E=210000;
Iz=4762;
L=1000;


% Give nodal displacements v1, theta1, v2, theta2: 
uu=[0 0 0 0]';      % TIP: Now all of nodal displacements are zero. You need to change this.      

% Let's visualize displacement field
x=0:L/100:L;

for ii=1:length(x)
    N=Shapef_EB4DOF(x(ii),L);
    uh(ii)=N*uu;
end

% Approximate displacement uh at x={500,750} with given u 
uh500=Shapef_EB4DOF(500,L)*uu
uh750=Shapef_EB4DOF(750,L)*uu

% Plots dislacement field 
figure(1)
plot(x,uh,'r-')
hold on;
plot(750,uh750,'b*',500,uh500,'b*')
xlabel('Longitudinal coordinate x [mm]')
ylabel('Displacment [mm]')

% Always check your internal forces that it does not produce any forces
% with zero displacements.
Fint=FintEB4DOF(E,Iz,L,uu)
