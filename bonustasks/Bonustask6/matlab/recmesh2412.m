function [P,nloc] = recmesh2412(n,m,a1,a2,b1,b2)
% makes a rectangular mesh with four node plate elements Q4 (2422)  

% generate nodal coordinates in xy-plane
% dx = a/n;
% dy = b/m;

% MKM 19082020: Modified to make it useful for visualization


dx = abs((a2-a1))/n;
dy = abs((b2-b1))/m;

% Geospace for spacing
xk=geospace(a1,a2,n+1,1)';
yk=geospace(b1,b2,m+1,1)';

P = [];
for k=1:m+1
    Pk = [xk yk(k).*ones(n+1,1)];  
    P = [P; Pk]; 
end

% generate elememnt connectivity
nloc = [];
for l = 1:1:m
  for k = 1:1:n
    loc = [k k+1 k+2+n n+k+1] + (l-1)*(n+1);
    nloc = [nloc; loc];
  end
end

