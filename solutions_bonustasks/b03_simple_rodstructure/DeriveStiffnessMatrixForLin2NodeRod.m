% BK10A6400 Basics of FE Analysis (FEMBasics2023)
% Teacher in charge: Marko Matikainen
% Code template for the Bonus Task 3.
%
% The code derives stiffness matrix for a two node linear rod (bar, truss) element 
% Coded by MKM for student's usage in the FEMBasics2023 course


%  Two node linear rod element ########################################
%   |-> x 
%   |<---- L ---->|  
%  
%  N1             N2 
%   o-------------o
%   ->            ->
%   u1            u2
%
%#######################################################################


% Lets use MATLAB's symbolic toolbox
clear all;
clc;
close all;


% Let's initiliaze variables
syms x L u1 u2 A E

uu=[u1 u2].';

AA=[1 0;
   1 L];

% polynomials
%p=[1,x].';               % Note that Matlab syms wants 
                         % notation .' for transpose. Also command
                         % transpose() works
p=transpose([1,x]);
N=p.'*AA^-1;

matlabFunction(N,'file','Shapef2NodeRod','vars',{x,L});

uh=N*uu;

% Axial strain
Epsxx=diff(uh,x);

% Internal strain energy for a rod element
Wintdx=1/2*E*A*Epsxx^2;

Wint=int(Wintdx,x,0,L)

% Fint= d Wint / d u
DOFs=2;

for kk=1:DOFs
    Fint(kk)=diff(Wint,uu(kk));
end

for ii=1:DOFs
    for kk=1:DOFs
        Kloc(ii,kk)=diff(Fint(ii),uu(kk));
    end
end

matlabFunction(Kloc,'file','KlocRod','vars',{A,E,L});






