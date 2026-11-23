function fe = EquivalentForcesBodyQ4(X,b,Lz)
% BK10A6400 Basics of FE Analysis (FEMBasics2026autumn)
% Teacher in charge: Marko Matikainen (MKM)
%
% Equivalent nodal forces from a constant body force for a general
% (distorted) Q4 element:
%
%   fe = int_V N^T b dV = Lz * int_{-1}^{1} int_{-1}^{1} N^T b det(Je) dxi deta
%
% For a Q4 element, det(Je) is linear in xi and eta, so the integrand is
% at most quadratic in each direction -> 2x2 Gauss quadrature is exact for
% a constant body force.
%
% Input:
%   X   : 8x1 nodal coordinates [X1;Y1;X2;Y2;X3;Y3;X4;Y4] (counterclockwise)
%   b   : 2x1 body force vector [N/m^3], e.g. b = rho*[0;-g]
%   Lz  : thickness [m]
% Output:
%   fe  : 8x1 equivalent nodal force vector [fx1;fy1;fx2;fy2;...]

Xn=X(1:2:end);  Yn=X(2:2:end);           % nodal x- and y-coordinates

gp=[-1 1]/sqrt(3);  w=[1 1];             % 2-point Gauss rule

fe=zeros(8,1);
for i=1:2
    for j=1:2
        xi=gp(i); eta=gp(j);
        [Nv,dNdxi,dNdeta]=ShapeQ4(xi,eta);
        % Jacobian matrix Je = [X,xi X,eta; Y,xi Y,eta]
        Je=[dNdxi*Xn, dNdeta*Xn;
            dNdxi*Yn, dNdeta*Yn];
        detJe=det(Je);
        if detJe<=0
            error('det(Je) <= 0: check node ordering (counterclockwise) or element shape.');
        end
        N=zeros(2,8);
        N(1,1:2:end)=Nv;  N(2,2:2:end)=Nv;
        fe=fe+w(i)*w(j)*Lz*(N.'*b)*detJe;
    end
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
