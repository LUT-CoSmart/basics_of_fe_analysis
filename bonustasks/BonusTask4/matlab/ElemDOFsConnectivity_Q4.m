function ElemDofsConnectivity= ElemDOFsConnectivity_Q4(ElemNodeConnectivity,NumberOfElems)
% Makes a tabel of elements' nodal displacements (DOFs) connectivity from elements' node connectivity table  

DofsAtElem=8;
ElemDofsConnectivity = zeros(NumberOfElems,DofsAtElem);
for k = 1:NumberOfElems
  n1 = [ElemNodeConnectivity(k,1)*2-1:1:ElemNodeConnectivity(k,1)*2];
  n2 = [ElemNodeConnectivity(k,2)*2-1:1:ElemNodeConnectivity(k,2)*2];
  n3 = [ElemNodeConnectivity(k,3)*2-1:1:ElemNodeConnectivity(k,3)*2];
  n4 = [ElemNodeConnectivity(k,4)*2-1:1:ElemNodeConnectivity(k,4)*2];
  ElemDofsConnectivity(k,:) = [n1 n2 n3 n4];
end
