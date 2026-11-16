function [NodeID,NumberOfIDs] = FindGlobNodeID(X0,x)
% Coded by MKM at LUT, make it more general later.
% Returns global nodal ID with given coordinate in x-direction
% and number of IDs
NodeID=[];
tol=2*sqrt(eps);
jj=1;
for ii=1:size(X0,1)
    if (abs(X0(ii,1)-x) < tol)
    %if (abs(X0(ii,1)-x) < tol) && (abs(X0(ii,2)-y) < tol)
        NodeID(jj)=ii;
        jj=jj+1;
    end    
end
NumberOfIDs=length(NodeID);
