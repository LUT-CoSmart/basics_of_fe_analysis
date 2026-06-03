% BK10A6400 Basics of FE Analysis (FEMBasics2025autumn)
% Teacher in charge: Marko Matikainen (MKM)
% Derivation of the bilinear plane element with symbolic volumetric
% integration
% Code template for the Bonus Task 5.

% Element's local numbering, let's defined a local coordinate system {x,y}
% to be at middle of a element
%
%
%  N4 ______ N3 
%    |      |
%    |      |   Ly
%    |______|
%   N1      N2
%      Lx 
%
%
%

clear all
close all
clc

DOFs=8;         % element has 4 nodes 2 dofs per node (4x2=8)

syms x y Lx Ly Lz E nu u1 u2 u3 u4 v1 v2 v3 v4

% A vector of nodal displacements
uu=[u1, v1, u2, v2, u3, v3, u4, v4].';
% A vector of nodal displacements separated into u and v
u=[u1, u2, u3, u4].';
v=[v1, v2, v3, v4].';

p=[1,x,y,x*y];

AA=zeros(4,4);

% TIPS: add components for the AA...if components are written you need 
% to check those!!!
AA=[1, 1/2*Lx, 1/2*Ly, 1/4*Lx*Ly;
    1, Lx, Ly, Lx*Ly;
    1, 1/4*Lx, -1/2*Ly, 1/4*Lx*Ly;
    1, -Lx, 1/2*Ly, 1/4*Lx*Ly]; % Add AA

% Shape functions
Nvec=p*AA^(-1);

% Writes shapefunctions into a file
matlabFunction(Nvec,'file','ShapefunctionsQ4','vars',{Lx,Ly,x,y});

% Displacement field interpolation in x and y directions
uh=Nvec*u;
vh=Nvec*v;

% TIPS: Define strains epsxx, epsyy and gammaxy...
% if they are (somehow) written, you need to check and fix those
% based on lectures and slides!!!
epsxx=diff(uh,x); % Eps_xx=partial uh / partial x
epsyy=diff(vh,y);
gammaxy=diff(uh,y)+diff(vh,x);

% Strain vector 
eps=[epsxx, epsyy, gammaxy].';

% Plain stress assumption
DD=E/(1-nu^2)*[1, nu, 0;
               nu, 1, 0;
               0, 0, (1-nu)/2];

% strain energy...note integration limits. Now the elemental 
% (local) coordinate system is located at x=0,y=0 (middle of element)
UdA=1/2*Lz*eps.'*DD*eps;
Udy=int(UdA,x,-Lx/2,Lx/2);          % integration over x analytically
U=int(Udy,y,-Ly/2,Ly/2);            % integration over y analytically


% Fe (size of 8)
% Internal forces, elastic forces etc.
for kk=1:DOFs,
     Fe(kk)=diff(U,uu(kk));
end       

% Stiffness matrix 8x8
for ii=1:DOFs,
    for jj=1:DOFs, 
        Kloc(ii,jj)=diff(Fe(ii),uu(jj));
    end    
end   


% Writes function for the elemental internal energy (strain energy). 
matlabFunction(U,'file','UQ4','vars',{E,nu,Lx,Ly,Lz,uu});

% Writes function for the elemental stiffness matrix. 
matlabFunction(Kloc,'file','KlocQ4','vars',{E,nu,Lx,Ly,Lz});