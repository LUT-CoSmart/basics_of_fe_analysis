% BK10A6400 Basics of FE Analysis (FEMBasics2024)
% Teacher in charge: Marko Matikainen 
% The code derives shape functions for a two node linear rod (bar, truss) element 
% Coded by MKM for student's usage in the FEMBasics2024 course

%   |<---- L ---->|  
%  
%  N1             N2 
%   o-------------o
%   ->            ->
%   u1            u2


clear all;
clc;
close all;

% Lets initiliaze variables
syms x L   
syms u1 u2 u3
% unknowns a_i
syms a0 a1 a2

% Vector of nodal coordinates
uu=[u1,u2,u3].';



% Polynomials
p=[1,x,x^2]';

A=[1 0 0;
   1 L/2 L^2/4;
   1 L L^2];

% Shape functions
N=zeros(3,1);
N=p'*A^-1;

% Writes function in a file
matlabFunction(N,'file','Shapef3NodeRod','vars',{x,L});




