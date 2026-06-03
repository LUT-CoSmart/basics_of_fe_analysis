function ID = DofID(DofsAtNode,nodes,comps)
%---------------------------------------------------------
% Postgraduate Course in Finite Element Method (FEM2019)
% Teacher in charge: Marko Matikainen
%---------------------------------------------------------
% xlocChosen.m  -  returns ID of nodal coordinates
% nodes  - List of Node IDs
% comps  - List of DOF's at node
% Coded by MKM for student's usage in the FEM2019 course
%---------------------------------------------------------
        
ID = [];
for n=1:length(nodes)
  nn = nodes(n);
  ID = [ID DofsAtNode*nn-DofsAtNode+comps];
end
  