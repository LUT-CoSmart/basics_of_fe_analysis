% BK10A6400 Basics of FE Analysis (FEMBasics2024autumn)
% Teacher in charge: Marko Matikainen (MKM)
% Template for the Bonus task 7 : Modal analysis. The code solves eigenvalues and modes by using Q4 elements.
%
% Goal: Derive lumped mass matrix, solve eigenfrequencies and make comparisoion with a consistent mass matrix. Solve the eigenfrequencies and modes of free, clamped and simply supported structures with two element types (KlocQ4xiFullGaussi and KlocQ4psR)
% Modal analysis for structure with mesh 1x1. Find shear and bending mode.
% Why eigenfrequencies are so different?
% Solve eight lowest eigenfrequencies and modes with mesh 


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


Analysis=1;
% number of modes visualized. If you use one element, maximum is 8
nmodesplotmax=8;

% Units are m,N
% Beam 1
E=1e9;
L=2;
H=0.4;
W=0.4;
nu=0.3;
rho=7850;


% Force (in vertical (global) direction at free end)
F=-62500;

%  Parameters for element and mesh
DofsAtNode=2;                               % 2 dofs at every node 
LocalSize=8;                                % Element's size (4 nodes * 2 Dofs at node)

% Mesh
% ScaleMesh=8;                                           % This you need to change. 
% nxElems=5*ScaleMesh;                                   % Discretization in x-direction      
% nyElems=1*ScaleMesh;                                 % Discretization in y-direction

ScaleMesh=1;                                           % This you need to change. 
nxElems=1*ScaleMesh;                                   % Discretization in x-direction      
nyElems=1*ScaleMesh;                                 % Discretization in y-direction



NumberOfElems=nxElems*nyElems;        
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

% Use DofID function to find out dof IDs for constraints
% DofId(DofsAtNode,[Nodes],[Dofs])
%bcDofIDs = DofID(2,[12,22],[1,2]); % Example about simply supported for 10x2 mesh
%bc(bcDofIDs)=0;

DOF=sum(bc);                                           % Number of unconstrained DOFs

AppliedForceNodeIDs=FindGlobNodeID(X0all,L);        % finds node IDs at free end (applied force) 
AppliedForceDofIDs=DofID(DofsAtNode,AppliedForceNodeIDs,[2]); 

% Initialization
Kglob=zeros(NumberOfAllDofs,NumberOfAllDofs);
Mglob=zeros(NumberOfAllDofs,NumberOfAllDofs);

% reshape makes vector form from X0all that is needed to create X in assembling 
X0allvec=reshape(X0all',[],1);
   
% Assemble the global matrix. This works for all cases
for ii=1:NumberOfElems
    X=X0allvec(ElemDOFsConnectivity(ii,:));
    

    % full Gaussian int
    %Kloc=KlocQ4xiFullGauss(X,E,nu,Lz);              
    % Selective reduced int
    Kloc=KlocQ4psR(X,E,nu,Lz);              % Function of nodal coordinates (trasnformationm via isoparametric)    
    %Mloc=MlocQ4(X,rho,Lz);                  % Consistent mass matrix
    Mloc=MlocQ4lumped(X,rho,Lz);           % Lumped


    for jj = 1:LocalSize   
        ind01 = ElemDOFsConnectivity(ii,jj); %Index 01
        for kk = 1:LocalSize
            ind02 = ElemDOFsConnectivity(ii,kk); % Index 02
            Kglob(ind01,ind02) = Kglob(ind01,ind02)+Kloc(jj,kk);
            Mglob(ind01,ind02) = Mglob(ind01,ind02)+Mloc(jj,kk); 
        end          
    end %End of Assembly
end


uglob=zeros(NumberOfAllDofs,1);       
fglob=zeros(NumberOfAllDofs,1);
fglob(AppliedForceDofIDs)=F/length(AppliedForceDofIDs);

% Boundary conditions (in a global coordinate systems)
Kglobred=zeros(DOF,DOF);
Mglobred=zeros(DOF,DOF);
fglobred=zeros(DOF,1);

% Boundary conditions. Eliminate linear constraints (by using bc vector)
Kglobred=Kglob(bc,bc);
Mglobred=Mglob(bc,bc);
fglobred=fglob(bc);



if Analysis==1,     % Eigenvalue
     disp('Modal analysis');
      
     [Eigenmodes,d] = eigs(Kglobred,Mglobred,nmodesplotmax,'smallestabs');
     Eigenfrequencies= real(sqrt(diag(d)))   % Eigenfrequencies of the gen. problem

     % Let's update modes with fixed coordinates (to visualize modes)
     EigenmodesFull=zeros(NumberOfAllDofs,nmodesplotmax);
     EigenmodesFull(bc,:)=Eigenmodes(:,1:nmodesplotmax);

  
elseif Analysis==2,          % Static
    disp('Static analysis'); 
    % Let's solve displacements uglob 
    uglobred=Kglobred\fglobred;
    
    % % Let's gather all displacements (solved and fixed)
    ugloball=zeros(NumberOfAllDofs,1);
    ugloball(bc)=uglobred;
else
    disp('****** Choose the analysis! ******');
end   


if Analysis==1,     % Eigenvalue
     % Post processing
    DrawingLimits=[-0.5 2.5 -0.5 1.5 -0.1 0.1]; % help axis
    
    figure(1)
    DrawMeshQ4(X0all,ElemNodeConnectivity,uglob,DofsAtNode,DrawingLimits,'b')             % initial conf
    % Node numbering
    % Writes nodal indeces in the figure
    for ii=1:NumberOfNodes,
        text(X0all(ii,1),X0all(ii,2),0, int2str(ii));            
    end    

    %DrawingLimits=[-150 50 -50 100 -10 10]; % help axis
    ScaleDisp=2;       % just to scale modes.
    
    % Undeformed and deformed
    for imode=1:nmodesplotmax,
        figure(imode+1)
        DrawMeshQ4(X0all,ElemNodeConnectivity,uglob,DofsAtNode,DrawingLimits,'b')             % initial conf
        DrawMeshQ4(X0all,ElemNodeConnectivity,EigenmodesFull(:,imode)*ScaleDisp,DofsAtNode,DrawingLimits,'r')   % deformed conf
        title(['Eigenfrequency: ' ,num2str(Eigenfrequencies(imode),5),' rad/s'],'FontName','Times New Roman','FontSize',[20]); 

        
        filenamefig=['ModalAnalysis_Mode',int2str(imode)];
        print('-dtiff', filenamefig);
        print('-depsc2', filenamefig);
        
    end
elseif Analysis==2,          % Static
    % Post processing
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

 
else
    disp('****** Choose the analysis! ******');
end   
