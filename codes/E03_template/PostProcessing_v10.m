% Postprocessing ********************************************************

% Displacements at applied loaded end
uFreeEnd=ugloball(AppliedForceDofIDs);
uFreeEndAverage= mean(uFreeEnd)

DrawingLimits=[80 220 30 120 -10 10];   % help axis
ScaleDisp=100;                          % just to scale displacement for visualization. If want to see realistic, use 1.

%Undeformed with node numbering. 
% TIP: From Fig, you can see what coordinate you need to use for function:
% FindGlobNodeID.
figure(1)
DrawMeshQ4(X0all,ElemNodeConnectivity,uglob,DofsAtNode,DrawingLimits,'b')             % initial conf
% Node numbering
%Writes nodal indeces in the figure
for ii=1:NumberOfNodes,
   text(X0all(ii,1),X0all(ii,2),0, int2str(ii));            
end


% Undeformed and deformed
figure(2)
DrawMeshQ4(X0all,ElemNodeConnectivity,uglob,DofsAtNode,DrawingLimits,'b')             % initial conf
DrawMeshQ4(X0all,ElemNodeConnectivity,ugloball*ScaleDisp,DofsAtNode,DrawingLimits,'r')   % deformed conf

    
% ***********************************************************************
% Next is not working properly in this code, do not use.
%PostProcessStresses_2412(X0all,ugloball,ElemDOFsConnectivity,ElemNodeConnectivity,NumberOfElems,ElemNodes,DofsAtNode,L,H,W,E,nu);
      