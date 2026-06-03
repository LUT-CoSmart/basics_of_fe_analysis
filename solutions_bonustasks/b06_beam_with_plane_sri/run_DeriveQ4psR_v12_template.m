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
%Xinit=[0,Lx,Lx,0].';
%Yinit=[0,0,Ly,Ly].';

Xinit=[X1,X3,X5,X7].';
Yinit=[X2,X4,X6,X8].';

Xh=Nvec*Xinit;
Yh=Nvec*Yinit;

% Strains
%epsxx=diff(uh,x);
%epsyy=diff(vh,y);
%gammaxy=diff(uh,y)+diff(vh,x);

% TIPS: define strains with binormalized coordinates xi. Check next ones does they
% make sense.
% ########################################################################
epsxx=diff(uh,xi)*diff(Xh,xi)^-1;
epsyy=diff(vh,eta)*diff(Yh,eta)^-1;
gammaxy=diff(uh,eta)*diff(Yh,eta)^-1+diff(vh,xi)*diff(Xh,xi)^-1;
%########################################################################

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
detJe=det(Je); 


DDs=DD(3,3);
DDn=DD(1:2,1:2);

% Lets split material stiffness matrix into shear and normal parts
DDs=DD(3,3);
DDn=DD(1:2,1:2);

UndA=1/2*eps(1:2).'*DDn*eps(1:2)*Lz*detJe;
UsdA=0;   % TIP: write shear part of energy. 


% Gauss integration
nip2=2;
[xiv2,wxi2]=gauleg2(-1,1,nip2);
[etav2,weta2]=gauleg2(-1,1,nip2);

Un=0;
for ii1=1:nip2
    for jj1=1:nip2
        Un=Un+subs(UndA,[xi,eta],[xiv2(ii1),etav2(jj1)])*wxi2(ii1)*weta2(jj1);
    end
end

nip1=1;   % TIP: set up a number of integartion points
[xiv,wxi]=gauleg2(-1,1,nip1);
[etav,weta]=gauleg2(-1,1,nip1);

Us=0;
for ii1=1:nip1
    for jj1=1:nip1
        Us=Us+subs(UsdA,[xi,eta],[xiv2(ii1),etav2(jj1)])*wxi2(ii1)*weta2(jj1);               % TIP: Write integration for shear energy. It should not be zero.
    end
end

U=Un+Us;


% strain energy...note integartion limits. Now the elemental 
% % (local) coordinate system is located at x=0,y=0 (middle of element)
% UdA=1/2*Lz*eps.'*DD*eps*det(Je);
% Udy=int(UdA,xi,-1,1);
% U=int(Udy,eta,-1,1);

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

%matlabFunction(Kloc,'file','KlocQ4psR','vars',{E,nu,Lx,Ly,Lz});
matlabFunction(Kloc,'file','KlocQ4psR','vars',{X,E,nu,Lz});         % no need for Lx and Ly anymore becuase they are ccomputed from X
