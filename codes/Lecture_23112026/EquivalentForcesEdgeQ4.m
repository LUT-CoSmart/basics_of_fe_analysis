function fe = EquivalentForcesEdgeQ4(X,edge,q,loadtype)
% BK10A6400 Basics of FE Analysis (FEMBasics2026autumn)
% Teacher in charge: Marko Matikainen (MKM)
%
% Equivalent nodal forces from a line load on one edge of a general
% (distorted) Q4 element:
%
%   fe = int_Gamma N^T q dGamma = int_{-1}^{1} N^T(s) q(s) J_edge(s) ds
%
% Edges are traversed counterclockwise with parameter s in [-1,1]:
%   edge 1: node 1 -> 2  (eta = -1, xi =  s)
%   edge 2: node 2 -> 3  (xi  =  1, eta = s)
%   edge 3: node 3 -> 4  (eta =  1, xi = -s)
%   edge 4: node 4 -> 1  (xi  = -1, eta = -s)
%
% Tangent vector:  tau = dX/ds = Je*c,   J_edge = norm(tau)
% Outward unit normal (counterclockwise node order): n = [tau_y; -tau_x]/J_edge
%
% Input:
%   X        : 8x1 nodal coordinates [X1;Y1;X2;Y2;X3;Y3;X4;Y4] (counterclockwise)
%   edge     : edge number 1..4
%   q        : load per unit length [N/m] (q = Lz*t)
%              'global'  : 2x1 constant [qx;qy], or 2x2 [q_start q_end]
%                          (linear variation from start node to end node of the edge)
%              'pressure': scalar p, or 1x2 [p_start p_end]; positive p
%                          pushes into the element (q = -p*n)
%   loadtype : 'global' (default) or 'pressure'
% Output:
%   fe       : 8x1 equivalent nodal force vector [fx1;fy1;fx2;fy2;...]
%
% 2-point Gauss rule is exact for a linearly varying load on a straight
% Q4 edge (N linear, q linear, J_edge constant).

if nargin<4, loadtype='global'; end

Xn=X(1:2:end);  Yn=X(2:2:end);

% Natural coordinates along the edge: [xi eta] = a + c*s
switch edge
    case 1, a=[0 -1]; c=[ 1  0];
    case 2, a=[1  0]; c=[ 0  1];
    case 3, a=[0  1]; c=[-1  0];
    case 4, a=[-1 0]; c=[ 0 -1];
    otherwise, error('edge must be 1, 2, 3 or 4');
end

% Load values at the start (s=-1) and end (s=1) of the edge
if strcmpi(loadtype,'pressure')
    if isscalar(q), q=[q q]; end
else
    if size(q,2)==1, q=[q q]; end
end

gp=[-1 1]/sqrt(3);  w=[1 1];             % 2-point Gauss rule

fe=zeros(8,1);
for k=1:2
    s=gp(k);
    xi=a(1)+c(1)*s;  eta=a(2)+c(2)*s;
    [Nv,dNdxi,dNdeta]=ShapeQ4(xi,eta);
    Je=[dNdxi*Xn, dNdeta*Xn;
        dNdxi*Yn, dNdeta*Yn];
    tau=Je*c.';                          % tangent dX/ds
    Jedge=norm(tau);                     % boundary Jacobian
    Ns=[(1-s)/2, (1+s)/2];               % linear interpolation of the load along the edge
    if strcmpi(loadtype,'pressure')
        n=[tau(2); -tau(1)]/Jedge;       % outward unit normal
        qs=-(Ns*q(:))*n;                 % pressure acts against the outward normal
    else
        qs=q*Ns.';                       % global components [qx;qy]
    end
    N=zeros(2,8);
    N(1,1:2:end)=Nv;  N(2,2:2:end)=Nv;
    fe=fe+w(k)*(N.'*qs)*Jedge;
end
end

function [Nv,dNdxi,dNdeta]=ShapeQ4(xi,eta)
% Bilinear shape functions, nodes 1(-1,-1) 2(1,-1) 3(1,1) 4(-1,1)
xin =[-1  1 1 -1];
etan=[-1 -1 1  1];
Nv    =0.25*(1+xin*xi).*(1+etan*eta);
dNdxi =0.25*xin.*(1+etan*eta);
dNdeta=0.25*(1+xin*xi).*etan;
end
