% Main code 
clear all; close all; clc;

% WRITE parameters E, A, L
E=1;A=1;L=1;

% Applied nodal force
Fx=1;

% WRITE a vector of applied exteral forces Fext
Fext=zeros(3,1);
Fext(3)=Fx;


% Call function
K1=Kloc(E,A,L/2)
K2=Kloc(E,A,L/2)

% Assembling
Kglob=zeros(3,3)

Kglob(1,1)=K1(1,1);
Kglob(1,2)=K1(1,2);
Kglob(2,1)=K1(2,1);
Kglob(2,2)=K1(2,2)+K2(1,1);

Kglob(3,2)=K2(2,1)
Kglob(3,3)=K2(2,2)
Kglob(2,3)=K2(1,2)
    

% Eliminate boundary conditions (bc)
%K1red=K1(2,2);
Kglobred=Kglob(2:3,2:3);
Fglobred=Fext(2:3);

% Solve displacement
u=Kglobred^-1*Fglobred

% Comparison with an analytical solution 
uanal=Fx*L/(E*A)


% Refresh your memory of basics constitutive eqs for one dimensional elongation problem (SLIDES)
 