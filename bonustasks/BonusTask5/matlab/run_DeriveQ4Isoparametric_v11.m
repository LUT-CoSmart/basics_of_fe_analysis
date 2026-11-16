% BK10A6400 Basics of FE Analysis (FEMBasics2024)
% Teacher in charge: Marko Matikainen (MKM)
% Derivation of the bilinear plane element
% Code template for the Bonus Task 6.

clear all
close all
clc

DOFs=8;         % element has 4 nodes 2 dofs per node (4x2=8)

syms xi eta Lx Ly Lz E nu 
syms u1 u2 u3 u4 v1 v2 v3 v4
syms X1 X2 X3 X4 X5 X6 X7 X8

% A vector of nodal displacements
uu=[u1, v1, u2, v2, u3, v3, u4, v4].';

% A vector of nodal position
X=[X1, X2, X3, X4, X5, X6, X7, X8].';

% A vector of nodal displacements separated into u and v
u=[u1, u2, u3, u4].';
v=[v1, v2, v3, v4].';

p=[1,xi,eta,xi*eta];

AA=zeros(4,4);

AA=[1, -1, -1, 1;
    1, 1, -1, -1;
    1, 1, 1, 1;
    1, -1, 1, -1];

Nvec=p*AA^(-1);

matlabFunction(Nvec,'file','ShapefunctionsQ4xi','vars',{Lx,Ly,xi,eta});

% Displacement field interpolation in x and y directions
uh=Nvec*u;
vh=Nvec*v;

% Interpolation for geometry (position)
% Xinit=[0,Lx,Lx,0].';
% Yinit=[0,0,Ly,Ly].';

Xinit=[X1,X3,X5,X7].';
Yinit=[X2,X4,X6,X8].';


Xh=Nvec*Xinit;
Yh=Nvec*Yinit;

% Strains
%epsxx=diff(uh,x);
%epsyy=diff(vh,y);
%gammaxy=diff(uh,y)+diff(vh,x);

epsxx=diff(uh,xi)*diff(Xh,xi)^-1;
epsyy=diff(vh,eta)*diff(Yh,eta)^-1;
gammaxy=diff(uh,eta)*diff(Yh,eta)^-1+diff(vh,xi)*diff(Xh,xi)^-1;


% Strain vector
eps=[epsxx, epsyy, gammaxy].';

% Plain stress
DD=E/(1-nu^2)*[1, nu, 0;
               nu, 1, 0;
               0, 0, (1-nu)/2];
Je(1,1)=diff(Xh,xi);
Je(1,2)=diff(Xh,eta);
Je(2,1)=diff(Yh,xi);
Je(2,2)=diff(Yh,eta);
           
           
% strain energy...note integartion limits. Now the elemental 
% (local) coordinate system is located at x=0,y=0 (middle of element)
UdA=1/2*Lz*eps.'*DD*eps*det(Je);
%Udy=int(UdA,xi,-1,1);
%U=int(Udy,eta,-1,1);

% Let's use numerical integrattion because derivation
% became very slow when X is included in derivation
% Gauss integration
nip2=2;
[xiv2,wxi2]=gauleg2(-1,1,nip2);
[etav2,weta2]=gauleg2(-1,1,nip2);

U=0;
for ii1=1:nip2
    for jj1=1:nip2
        U=U+subs(UdA,[xi,eta],[xiv2(ii1),etav2(jj1)])*wxi2(ii1)*weta2(jj1);
    end
end


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

%matlabFunction(Kloc,'file','KlocQ4xi','vars',{E,nu,Lx,Ly,Lz});
matlabFunction(Kloc,'file','KlocQ4xi','vars',{X,E,nu,Lz});         % no need for Lx and Ly anymore becuase they are ccomputed from X

