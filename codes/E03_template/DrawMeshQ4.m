function h = DrawMeshQ4(X0,nloc,u,DofsAtNode,DrawingLimits,color)
% draws mesh from node to node.
set(0,'DefaultLineLineWidth',1.5)

[NumberOfNodes,~] = size(X0);
[NumberOfElems,~] = size(nloc);
NumberOfDofs = DofsAtNode*NumberOfNodes;

% Makes an array [NumberOfNodes,NumberOfTransDofs]
% for drawing from displacement vector
loc = DofID(DofsAtNode,1:NumberOfNodes,[1 2]); % choose trans dofs
Nu = reshape(u(loc),DofsAtNode,NumberOfNodes).';
C = X0 + Nu;     % updates X0 with nodal displacments

% Instructions for drawing lines from node to node (local numbering)
nlocdraw=[1,2,3,4,1];

for ii=1:NumberOfElems
    Cx=C(nloc(ii,nlocdraw),1);
    Cy=C(nloc(ii,nlocdraw),2);
    Cz=zeros(5,1);
    plot3(Cx,Cy,Cz,color)
    hold on
end

axis equal
axis(DrawingLimits);
view(2)