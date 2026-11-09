% Main code 
clear all; close all; clc;

% WRITE parameters E, A, L
E=210000;A=pi*10^2;L=3000;

% Applied nodal force
Fx=1000;

% WRITE a vector of applied exteral forces Fext
Fext=zeros(2,1);
Fext(2)=Fx;


% Call function
K1=Kloc(E,A,L)
Kglob=K1;


    

% Eliminate boundary conditions (bc)
Kglobred=K1(2,2);
Fextglobred=Fext(2);


% Solve displacement
u=Kglobred^-1*Fextglobred

% Comparison with an analytical solution 
uanal=Fx*L/(E*A)


% Refresh your memory of basics constitutive eqs for one dimensional elongation problem (SLIDES)
 