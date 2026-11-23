% BK10A6400 Basics of FE Analysis (FEMBasics2026autumn)
% Teacher in charge: Marko Matikainen (MKM)
% Derivation of a vector of elemental equivalent forces for Q4 element
%
% NOTE: The element geometry is a rectangle Lx x Ly aligned with the global
% x- and y-axes (Xinit, Yinit below). The generated functions are therefore
% valid ONLY for such undistorted rectangular elements. For a general
% (distorted) Q4 element, detJe must be evaluated from the actual nodal
% coordinates and the line element along an edge is
%   ds = sqrt((dX/dxi)^2 + (dY/dxi)^2) dxi   (edges eta = const.)
%   ds = sqrt((dX/deta)^2 + (dY/deta)^2) deta (edges xi = const.)
%
%   Element definition (natural coordinates xi, eta in [-1,1])
%
%          N4 (-1,1)    N3 (1,1)
%             o-----------o
%             |           |
%             |           |
%             o-----------o
%          N1 (-1,-1)   N2 (1,-1)
%
% Nodal dof order: [u1, v1, u2, v2, u3, v3, u4, v4]

clear all
close all
clc

syms xi eta Lx Ly Lz

% Shape functions from the polynomial basis p = [1, xi, eta, xi*eta]
p=[1,xi,eta,xi*eta];

% Rows: basis p evaluated at nodes N1..N4
AA=[1, -1, -1,  1;
    1,  1, -1, -1;
    1,  1,  1,  1;
    1, -1,  1, -1];

Nvec=p*AA^(-1);
matlabFunction(Nvec,'file','ShapefunctionsQ4xi','vars',{Lx,Ly,xi,eta});

% Interpolation for geometry (position), rectangular element Lx x Ly
Xinit=[0,Lx,Lx,0].';
Yinit=[0,0,Ly,Ly].';

% Position field
Xh=Nvec*Xinit;
Yh=Nvec*Yinit;

% Jacobian matrix (here Je = [dX/dxi dX/deta; dY/dxi dY/deta],
% i.e. the transpose of the usual convention; the determinant is the same)
Je(1,1)=diff(Xh,xi);
Je(1,2)=diff(Xh,eta);
Je(2,1)=diff(Yh,xi);
Je(2,2)=diff(Yh,eta);
detJe=det(Je);          % = Lx*Ly/4 for the rectangle

% Shape function matrix (2x8) in xi coordinate system
% Row 1 -> u-dofs (columns 1,3,5,7), row 2 -> v-dofs (columns 2,4,6,8)
N=sym(zeros(2,8));
for ii=1:2
 for jj=1:4
     jj2=(ii-1)+(jj-1)*2+1;
     N(ii,jj2)=Nvec(jj);
 end
end

% Example 1: structure's own mass as loading.
% Equivalent force vector (element level) from the structure's own weight.
% Gravity acts in the negative y-direction.
syms g rho
gg=[0,-g].';     % a field of gravity
b=rho*gg;        % body force per unit volume
fequivBF_dxideta=Lz*b.'*N*detJe;
fequivBF=int(int(fequivBF_dxideta,xi,-1,1),eta,-1,1)
matlabFunction(fequivBF,'file','EquivalentForcesBF','vars',{Lx,Ly,Lz,g,rho});

% Example 2: Equivalent nodal forces from a constant line load on the top
% edge (eta = 1, nodes N3 and N4). The load q acts DOWNWARDS (-y).
% q = Lz*t, i.e. force per unit length (t = surface traction).
% Along the edge eta = 1: ds = dX/dxi dxi = Je(1,1) dxi
syms q
qq=[0,-q].';
fequivConstLineLoad_dxi=subs(qq.'*N*Je(1,1),eta,1);
fequivConstLineLoad=int(fequivConstLineLoad_dxi,xi,-1,1)
matlabFunction(fequivConstLineLoad,'file','EquivalentForcesConstantLL','vars',{Lx,q});

% Example 3: Equivalent nodal forces from a constant line load on the bottom
% edge (eta = -1, nodes N1 and N2). The load q acts UPWARDS (+y).
% q = Lz*t, i.e. force per unit length.
% Along the edge eta = -1: ds = dX/dxi dxi = Je(1,1) dxi
qq=[0,q].';
fequivConstLineLoadBottom_dxi=subs(qq.'*N*Je(1,1),eta,-1);
fequivConstLineLoadBottom=int(fequivConstLineLoadBottom_dxi,xi,-1,1)
matlabFunction(fequivConstLineLoadBottom,'file','EquivalentForcesConstantLLBottom','vars',{Lx,q});

% Example 4: Equivalent nodal forces from a constant line load on the right
% side (xi = 1, nodes N2 and N3) in the longitudinal direction.
% The load q acts in the -x direction (i.e. compression on the right edge).
% q = Lz*t, i.e. force per unit length.
% Along the edge xi = 1: ds = dY/deta deta = Je(2,2) deta
qq=[-q,0].';
fequivConstLineLoadSide_dxi=subs(qq.'*N*Je(2,2),xi,1);
fequivConstLineLoadSide=int(fequivConstLineLoadSide_dxi,eta,-1,1)
matlabFunction(fequivConstLineLoadSide,'file','EquivalentForcesConstantLLSide','vars',{Ly,q});

% Expected results (rectangular element):
%   Own weight : v1 = v2 = v3 = v4 = -rho*g*Lx*Ly*Lz/4
%   Top        : v3 = v4 = -q*Lx/2
%   Bottom     : v1 = v2 = +q*Lx/2
%   Side       : u2 = u3 = -q*Ly/2
% Note: the vectors are returned as 1x8 row vectors (b.'*N).
