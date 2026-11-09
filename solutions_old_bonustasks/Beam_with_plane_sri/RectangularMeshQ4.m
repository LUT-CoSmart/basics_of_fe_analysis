function [X0all,NodeConnect] = RectangularMeshQ4(n,m,a1,a2,b1,b2)
% makes a rectangular mesh with four node plate elements Q4  
% generate nodal coordinates in xy-plane
% dx = a/n;
% dy = b/m;

% MKM 21032022: Modified for the linear FE course

dx = abs((a2-a1))/n;
dy = abs((b2-b1))/m;

% Geospace for spacing
xk=geospace(a1,a2,n+1,1)';
yk=geospace(b1,b2,m+1,1)';

X0all = [];
for k=1:m+1
    X0allk = [xk yk(k).*ones(n+1,1)];  
    X0all = [X0all; X0allk]; 
end

% generate elememnt connectivity
NodeConnect = [];
for l = 1:1:m
  for k = 1:1:n
    NodeConnectk = [k k+1 k+2+n n+k+1] + (l-1)*(n+1);
    NodeConnect = [NodeConnect; NodeConnectk];
  end
end