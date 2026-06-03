% BK10A6400 Basics of FE Analysis (FEMBasics2025)
% Teacher in charge: Marko Matikainen (MKM)
% Template for the bonus task 6. This is quiet similar with a template for the Exercise 3. The code solves global displacements at free end 
% by using Q4 elements.
%
% Goal: Solve the problem with two elements (KlocQ4xi and KlocQ4psR)
% and visualize the rate of convergence using those two elements
% Coded by MKM for student's usage in the FEMBasics2025 course

%   Element definition 
%
%          N4     N3      
%          
%          \            
%          \o-----o
%          \|     |     
%          \|     |     
%          \o-----o
%          \
%         N1     N2      

% First automatized matlab-based FE code for a simple beam structure analyzed by Q4
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
LocalSize=8;                                % Element's size (4 nodes * 2 Dofs at node)

% Mesh
nxElems=4;                                  % This you need to change. Discretization in x-direction      
nyElems=1;                                  % Discretization in y-direction
NumberOfElems=nxElems*nyElems;

Lx=L/nxElems;                               % Element's size in x-direction
Ly=H/nyElems;                               % Element's size in y-direction                                 
Lz=W;                                 

%ElemNodeConnectivity=[1,2,4,3];             % for one element. You need two lines for two elements...
% Ele1 N1 N2 N5 N4  u1,v1,u2,v2,u5,v5,u4,v4

% Returns nodal position coordinates X0all in initial confoguration and element connectivity
[X0all,ElemNodeConnectivity]=RectangularMeshQ4(nxElems,nyElems,0,L,0,H);
NumberOfNodes=size(X0all,1);
NumberOfAllDofs = DofsAtNode*NumberOfNodes; % Number of all dofs of structure

% Lets use function that creates DOFs connectivity
ElemDOFsConnectivity=ElemDOFsConnectivity_Q4(ElemNodeConnectivity,NumberOfElems);

% Let's define bc vector (full of ones and zeros) to set up boundary
% conditions (now linear constraints only) easier
bc = logical(ones(1,NumberOfAllDofs));                  % 1 = no constraints, 0 = constraints
bcNodeIDs = FindGlobNodeID(X0all,0);                    % let's find global node ID's x=0
bcDofIDs=DofID(DofsAtNode,bcNodeIDs,[1:2]);             % fixed coordinates ID          
bc(bcDofIDs)=0;
DOF=sum(bc);                                           % Number of unconstrained DOFs

AppliedForceNodeIDs=FindGlobNodeID(X0all,L);        % finds node IDs at free end (applied force) 
AppliedForceDofIDs=DofID(DofsAtNode,AppliedForceNodeIDs,[2]); 

% Initialization
Kglob=zeros(NumberOfAllDofs,NumberOfAllDofs);

% reshape makes vector form from X0all that is needed to create X in assembling 
X0allvec=reshape(X0all',[],1);
   
% Assemble the global matrix. This works for all cases
for ii=1:NumberOfElems
    X=X0allvec(ElemDOFsConnectivity(ii,:));
    % full Gaussian int
    Kloc=KlocQ4xi(X,E,nu,Lz);              
    % Selective reduced int
    %Kloc=KlocQ4psR(X,E,nu,Lz);              % Function of nodal coordinates (trasnformationm via isoparametric)    
 
    for jj = 1:LocalSize   
        ind01 = ElemDOFsConnectivity(ii,jj); %Index 01
        for kk = 1:LocalSize
            ind02 = ElemDOFsConnectivity(ii,kk); % Index 02
            Kglob(ind01,ind02) = Kglob(ind01,ind02)+Kloc(jj,kk);           
        end          
    end %End of Assembly
end


uglob=zeros(NumberOfAllDofs);       
fglob=zeros(NumberOfAllDofs,1);
fglob(AppliedForceDofIDs)=F/length(AppliedForceDofIDs);

% Boundary conditions (in a global coordinate systems)
Kglobred=zeros(DOF,DOF);
fglobred=zeros(DOF,1);

% Boundary conditions. Eliminate linear constraints (by using bc vector)
Kglobred=Kglob(bc,bc);
fglobred=fglob(bc);

% Let's solve displacements uglob 
uglobred=Kglobred\fglobred;

% % Let's gather all displacements (solved and fixed)
ugloball=zeros(NumberOfAllDofs,1);
ugloball(bc)=uglobred;

% Let's present averaged displacement at free end
uglobVerticalFreeEnd=mean(ugloball(AppliedForceDofIDs)) 




DrawingLimits=[-0.5 2.5 -0.5 1.5 -0.1 0.1]; % help axis

%DrawingLimits=[-150 50 -50 100 -10 10]; % help axis
ScaleDisp=100;       % just to scale displacement. If want to see realistic, use 1.

% Undeformed and deformed
figure(1)
DrawMeshQ4(X0all,ElemNodeConnectivity,uglob,DofsAtNode,DrawingLimits,'b')             % initial conf
DrawMeshQ4(X0all,ElemNodeConnectivity,ugloball*ScaleDisp,DofsAtNode,DrawingLimits,'r')   % deformed conf

% Undeformed with node numbering
figure(2)
DrawMeshQ4(X0all,ElemNodeConnectivity,uglob,DofsAtNode,DrawingLimits,'b')             % initial conf
% Node numbering
% Writes nodal indeces in the figure
for ii=1:NumberOfNodes,
    text(X0all(ii,1),X0all(ii,2),0, int2str(ii));            
end    
% ***********************************************************************
