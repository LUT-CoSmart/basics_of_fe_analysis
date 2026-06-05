% BK10A6400 Basics of FE Analysis (FEMBasics2025autumn)
% Teacher in charge: Marko Matikainen (MKM)
% Template for the Bonus Task 5.
%
% Goal: The code solves global displacements at free end by using a bilinear four node element. 
% Coded by MKM for student's usage in the FEMBasics2024 course


%          N4     N5    N6    
%          
%          \            |F
%          \o-----o-----o
%          \|     |     |
%          \|     |     |
%          \o-----o-----o
%          \
%         N1     N2    N3  


%          N3     N4      
%          
%          \            
%          \o-----o
%          \|     |     
%          \|     |     
%          \o-----o
%          \
%         N1     N2      


% Ele1 N1 N2 N4 N3  
% u1,v1,u2,v2,u4,v4,u3,v3

% So-called connectivity matrix
% Ele1 N1 N2 N5 N4  u1,v1,u2,v2,u5,v5,u4,v4
% Ele2 N2 N3 N6 N5  u2,v2,u3,v3,u6,v6,u5,v5,
% so nodes 2 and 5 are common nodes

% First matlab-based FE code for a simple beam structure analyzed by Q4
clear all;
close all;
clc;
format long;

% Units are m,N
% Beam 1
E=2.07e11;
L=2;
H=0.5;
W=0.1;
nu=0.3;

% Force (in vertical (global) direction at free end)
F=-62500;

%  Parameters for element and mesh
DofsAtNode=2;                               % 2 dofs at every node 
nxElems=2;                                  % TIP: For 2*1 mesh, you need to change this. Discretization in x-direction      
nyElems=1;                                  % Discretization in y-direction
NumberOfElems=nxElems*nyElems;
%NumberOfNodes=4;
NumberOfNodes=6;                           % TIP: For 2*1 mesh, you need to uncomment this.    
NumberOfAllDofs = DofsAtNode*NumberOfNodes; % Number of all dofs of structure

Lx=L/nxElems;                               % Element's size in x-direction
Ly=H/nyElems;                               % Element's size in y-direction                                 
Lz=W;                                 

ElemNodeConnectivity=[1,2,5,4;
                    2 3 6 5];             % TIP: for one element. You need two lines for two elements...
                                            % TIP: so it goes something
                                            % like: ElemNodeConnectivity=[1,2,5,4; 
                                            %                             ?,?,?,?];  

% Ele1 N1 N2 N5 N4  u1,v1,u2,v2,u5,v5,u4,v4
% Ele2 N2 N3 N6 N5 

% Lets use functions that creates DOFs connectivity
ElemDOFsConnectivity=ElemDOFsConnectivity_Q4(ElemNodeConnectivity,NumberOfElems);

% Elemental stiffness matrix (all are now same size)
Kloc=KlocQ4(E,nu,Lx,Ly,Lz);

% Initialization
Kglob=zeros(NumberOfAllDofs,NumberOfAllDofs);


% Assemble the global matrix. This works for all cases
localsize=size(Kloc,1);
for ii=1:NumberOfElems 
for jj = 1:localsize    
    ind01 = ElemDOFsConnectivity(ii,jj); %Index 01
    for kk = 1:localsize
        ind02 = ElemDOFsConnectivity(ii,kk); % Index 02
        Kglob(ind01,ind02) = Kglob(ind01,ind02)+Kloc(jj,kk);           
    end          
 end %End of Assembly
end

% % One element
% % Load vector (in a global coordinate systems)
fglob=zeros(12,1);                                               % TIP: Initialize a vector of external forces. If system has 2x1 elements, it means 6 nodes and 2 nodal coordinates at a node...
% Note you need to add loading with respect to correct dofs
fglob(6)=F/2;        % add loading                              % TIP: Modify so that load is applied at the free end in case of 2x1 mesh. So force is applied in vertical directions at node 3.
fglob(12)=F/2;        % add loading                              % TIP: Modify so that load is applied at the free end in case of 2x1 mesh. So force is applied in vertical directions at node 6.

% Boundary conditions (in a global coordinate systems)
% Uncomment and add correct numbers 
Kglobred=zeros(8,8);                                            % TIP: Initialize a reduced global stiffness matrix (constraints eliminated). If all nodal displacements at support are fixed, then you have four constraints...
fglobred=zeros(8,1);                                            % TIP: Initialize a vector of external forces (constraints eliminated).
Kglobred=Kglob([3:6,9:12],[3:6,9:12]);                            % TIP: Eliminated rows and columns regarding constrained nodal displacements. Here I just picked up unconstrained columns and rows.                           
fglobred=fglob([3:6,9:12]);                                      % TIP: Eliminated rows regarding constrained nodal displacements. Here I just picked up unconstrained rows.   

% Let's solve displacements uglob 
uglobred=Kglobred\fglobred;    % uglobred=Kglobred^-1*fglobred;

% % Let's gather all displacements (solved and fixed)
ugloball=zeros(12,1);                                            % TIP: ugloball includes all nodal displacements (constrained (fixed) and free ones). Make initiliazation.
ugloball([3:6,9:12])=uglobred                                    % TIP: update ugloball with solved nodal displacements. 





