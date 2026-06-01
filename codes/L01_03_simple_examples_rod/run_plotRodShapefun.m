clear all;
close all;
clc;

L=2;

% nodal displacement vectors from the assigment
u0=[0,0]';
udef=[0.1,0.1]';

x=0:0.1:L;
for ii=1:length(x),
    N=Shapef2NodeRod(x(ii),L);    
    N1(ii)=N(1);
    N2(ii)=N(2);
    uh0(ii)=N*u0;       % displacement field with u=[0,0]'
    uhdef(ii)=N*udef;   % displacement field with u=[0,0.1]';
end


% Plots shape functions and save figure
figure(1)
plot(x,N1,'r-')
hold on;
plot(x,N2,'b-')
legend('N_1','N_2','Location','North')
xlabel('Longitudinal coordinate x [m]')
ylabel('Shapefunction')
print("Shapef2NodeRod",'-depsc2');
print("Shapef2NodeRod",'-dtiff');
set(gca,'FontSize',18,'FontName','Times');

% Plots displacemnet field
figure(2)
plot(x,uh0,'r-')
hold on;
plot(x,uhdef,'b-')
legend('Displacement at u=[0,0]','Displacement at u=[0,0.1]','Location','North')
xlabel('Longitudinal coordinate x [m]')
ylabel('Displacement [m]')
print("DisplacementField",'-depsc2');
print("DisplacementField",'-dtiff');
set(gca,'FontSize',18,'FontName','Times');


