function [NodalID,NumberOfIDs] = FindGlobNodalID(P,x)
% Coded by MKM, make it more general later.
% Returns global nodal ID with given coordinate in x-direction
% and number of IDs
NodalID=[];
tol=2*sqrt(eps);
jj=1;
for ii=1:size(P,1)
    if (abs(P(ii,1)-x) < tol)
    %if (abs(P(ii,1)-x) < tol) && (abs(P(ii,2)-y) < tol)
        NodalID(jj)=ii;
        jj=jj+1;
    end    
end
NumberOfIDs=length(NodalID);
