% BK10A6400 Basics of FE Analysis (FEMBasics2026)
% Teacher in charge: Marko Matikainen (MKM)
% Derivation of the 8-node quadratic (serendipity) plane element, Q8
% Code template for the Bonus Task 5.
clear all
close all
clc
DOFs=16;        % element has 8 nodes 2 dofs per node (8x2=16)
syms xi eta Lz E nu
syms u1 u2 u3 u4 u5 u6 u7 u8 v1 v2 v3 v4 v5 v6 v7 v8
syms X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15 X16
% Node numbering: corners counter-clockwise, then midside nodes
%   4 --- 7 --- 3
%   |           |
%   8           6
%   |           |
%   1 --- 5 --- 2
% A vector of nodal displacements
uu=[u1, v1, u2, v2, u3, v3, u4, v4, u5, v5, u6, v6, u7, v7, u8, v8].';
% A vector of nodal position [x1 y1 x2 y2 ... x8 y8]
X=[X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14, X15, X16].';
% A vector of nodal displacements separated into u and v
u=[u1, u2, u3, u4, u5, u6, u7, u8].';
v=[v1, v2, v3, v4, v5, v6, v7, v8].';
% Polynomial basis of the serendipity element (no xi^2*eta^2 term)
p=[1, xi, eta, xi^2, xi*eta, eta^2, xi^2*eta, xi*eta^2];
% Each row: p evaluated at one node (xi_i, eta_i)
%      1  xi  eta xi^2 xi*eta eta^2 xi^2*eta xi*eta^2
AA=[1, -1, -1,  1,   1,    1,    -1,      -1;    % node 1 (-1,-1)
    1,  1, -1,  1,  -1,    1,    -1,       1;    % node 2 ( 1,-1)
    1,  1,  1,  1,   1,    1,     1,       1;    % node 3 ( 1, 1)
    1, -1,  1,  1,  -1,    1,     1,      -1;    % node 4 (-1, 1)
    1,  0, -1,  0,   0,    1,     0,       0;    % node 5 ( 0,-1)
    1,  1,  0,  1,   0,    0,     0,       0;    % node 6 ( 1, 0)
    1,  0,  1,  0,   0,    1,     0,       0;    % node 7 ( 0, 1)
    1, -1,  0,  1,   0,    0,     0,       0];   % node 8 (-1, 0)
Nvec=simplify(p*AA^(-1));
matlabFunction(Nvec,'file','ShapefunctionsQ8xi','vars',{xi,eta});
% Check: shape functions sum to one (partition of unity)
assert(isequal(simplify(sum(Nvec)),sym(1)))
% Displacement field interpolation in x and y directions
uh=Nvec*u;
vh=Nvec*v;
% Interpolation for geometry (position), isoparametric
Xinit=X(1:2:end);     % [X1 X3 ... X15].'
Yinit=X(2:2:end);     % [X2 X4 ... X16].'
Xh=Nvec*Xinit;
Yh=Nvec*Yinit;
% Jacobian matrix  Je = d(x,y)/d(xi,eta)
Je=[diff(Xh,xi), diff(Xh,eta);
    diff(Yh,xi), diff(Yh,eta)];
detJe=det(Je);
% Strains. Chain rule:  [d/dxi; d/deta] = Je.' * [d/dx; d/dy]
%   =>  [d/dx; d/dy] = inv(Je.') * [d/dxi; d/deta]
% (valid for any element shape, also distorted / curved elements)
invJT=inv(Je.');
du=invJT*[diff(uh,xi); diff(uh,eta)];   % [du/dx; du/dy]
dv=invJT*[diff(vh,xi); diff(vh,eta)];   % [dv/dx; dv/dy]
epsxx=du(1);
epsyy=dv(2);
gammaxy=du(2)+dv(1);
% Strain vector
eps=[epsxx, epsyy, gammaxy].';
% Plane stress
DD=E/(1-nu^2)*[1, nu, 0;
               nu, 1, 0;
               0, 0, (1-nu)/2];
% Strain-displacement matrix B (eps = B*uu, eps is linear in uu)
BB=jacobian(eps,uu);                    % 3x16
% Strain energy U = 1/2*Lz*int(eps.'*DD*eps*detJe) and K = d^2U/duu^2.
% Because eps = BB*uu, this is  K = Lz*int(BB.'*DD*BB*detJe), which is
% much faster to evaluate than differentiating U twice symbolically.
% Gauss integration, full integration for Q8: 3x3 points
nip3=3;
[xiv3,wxi3]=gauleg2(-1,1,nip3);
[etav3,weta3]=gauleg2(-1,1,nip3);
Kloc=sym(zeros(DOFs));
for ii1=1:nip3
    for jj1=1:nip3
        Bg=subs(BB,[xi,eta],[xiv3(ii1),etav3(jj1)]);
        detJg=subs(detJe,[xi,eta],[xiv3(ii1),etav3(jj1)]);
        Kloc=Kloc+Lz*(Bg.'*DD*Bg)*detJg*wxi3(ii1)*weta3(jj1);
    end
end
% Internal (elastic) forces, if needed: Fe = Kloc*uu
% Stiffness matrix 16x16
matlabFunction(Kloc,'file','KlocQ8xi','vars',{X,E,nu,Lz});
% Quick check with a rectangular element (Lx = 2, Ly = 1):
% Xr=[0 0, 2 0, 2 1, 0 1, 1 0, 2 0.5, 1 1, 0 0.5].';
% K=KlocQ8xi(Xr,210e9,0.3,0.01);
% norm(K-K.')/norm(K)    % ~0  (symmetric)
% sum(abs(eig(K))<1e-6*max(abs(eig(K))))   % = 3 rigid-body modes