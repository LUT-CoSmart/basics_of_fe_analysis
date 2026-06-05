% BK10A6400 Basics of FE Analysis (FEMBasics2022)
% Teacher in charge: Marko Matikainen (MKM)
% Derivation of the bilinear plane element
% Code template for the Bonus Task 8.

clear all
close all
clc

DOFs=8;         % element has 4 nodes 2 dofs per node (4x2=8)

syms xi eta Lx Ly Lz E nu 
syms u1 u2 u3 u4 v1 v2 v3 v4
syms X1 X2 X3 X4 X5 X6 X7 X8
syms rho

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


Je(1,1)=diff(Xh,xi);
Je(1,2)=diff(Xh,eta);
Je(2,1)=diff(Yh,xi);
Je(2,2)=diff(Yh,eta);
detJe=det(Je); 


% Consistent mass matrix (via shape functions)
% Shape function matrix in x coordinate system
%Nm=zeros(8,8);
    for ii=1:2
     for jj=1:4
         jj2=(ii-1)+(jj-1)*2+1;
         N(ii,jj2)=Nvec(jj);
     end
    end
   
NNdxideta=N.'*N*detJe;
NNdeta=int(NNdxideta,xi,-1,1);
NN_intdxideta=int(NNdeta,eta,-1,1);
Mloc=rho*Lz*NN_intdxideta;

matlabFunction(Mloc,'file','MlocQ4','vars',{X,rho,Lz});       

% Lumped mass matrix
Mlumpedtmp=zeros(8,8);
alpha=1;
Mlumpedtmp=alpha*rho*Lz*diag(ones(8,1)); % TIP: you need to think how much mass should ...
%be lumped at nodes...so this line must be modified by multiplying it by the ...
%factor alpha. Should be different than 1.
detJeint=int(int(detJe,xi,-1,1),eta,-1,1);
Mloclumbed=Mlumpedtmp*detJeint;
matlabFunction(Mloclumbed,'file','MloclumpedQ4','vars',{X,rho,Lz});




% MdV=rho*N.'*N*Lz*detJe;
% 
% M=zeros(8,8);
% for ii=1:8,
%    for jj=1:8,
%       M(ii,jj)=int(int(MdV(ii,jj),xi,-1,1),eta,-1,1);
%    end
% end




%M=zeros(8,8);
%for ii=1:8,
%    for jj=1:8,
%        for ii1=1:nip2
%            for jj1=1:nip2
%                M(ii,jj)=M(ii,jj)+subs(MdV(ii,jj),[xi,eta],[xiv2(ii1),etav2(jj1)])*wxi2(ii1)*weta2(jj1);
%            end
%        end
%    end
%end




