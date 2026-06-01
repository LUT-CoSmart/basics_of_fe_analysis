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
syms u1 u2
% unknowns a_i
syms a0 a1

% Vector of nodal coordinates
uu=[u1,u2].';



% Polynomials
p=[1,x]';

A=[1 0;
   1 L];

% Shape functions
N=zeros(2,1);
N=p'*A^-1;

% Writes function in a file
matlabFunction(N,'file','Shapef2NodeRod','vars',{x,L});




