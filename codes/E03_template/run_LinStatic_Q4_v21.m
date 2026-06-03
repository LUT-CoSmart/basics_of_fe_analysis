% BK10A6400 Basics of FE Analysis (FEMBasics2023)
% Teacher in charge: Marko Matikainen (MKM)
% Template for the Exercise 3.
% Plate with a hole 
% Solves also simple cantilever probelm using three different meshes.
% Choose one you want to use.
% 
% Goal: The code solves global displacements at free end by using Q4 element. 
% Coded by MKM for student's usage in the FEMBasics2023 course
%
% You need to set up loading and boundary conditions by using double symmetry              
% Now loadings and bc's are not correctly applied.
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

% Automatized matlab-based linear FE code for a simple beam structure with hole analyzed by Q4
clear all;
close all;
clc;
format long;

% Units are mm,N,Mpa
% Parameters for a plate wth hole
% These parameters are for the double symmetric part.
E=210000;
L=100;
H=50;
W=5;
nu=0.3;

% Force (in vertical (global) direction at free end)
%F=25000;
q=100;  % [N/mm^2]
F=H*W*q;    % 5000;   % Total force 5000 N: q=100 N/mm2, F=q*H*W;

%  Parameters for element and mesh
DofsAtNode=2;                               % 2 dofs at every node 
ElemNodes=4;
LocalSize=ElemNodes*DofsAtNode;                                % Element's size (4 nodes * 2 Dofs at node)
Mesh=100;                               % Choose 100 or 800 elements.
% Element's size in y-direction                                 
Lz=W;                                 

% Read the matrix from the saved file using readmatrix

if Mesh==100,   % 100 elements
    NLIST= readmatrix('NLIST_mesh100.txt');      % node locations in x,y coordinates 
    ELIST= readmatrix('ELIST_mesh100.txt');      % element's connectivity matrix
elseif Mesh==800 % 800 elements
    NLIST= readmatrix('NLIST_mesh800.txt');      % node locations in x,y coordinates 
    ELIST= readmatrix('ELIST_mesh800.txt');      % element's connectivity matrix
else

end    

X0all(:,1:2)=NLIST(:,2:3);
ElemNodeConnectivity(:,1:4)=ELIST(:,[2,5,4,3]);      % FEMAP defines elements in opposiet directions?!
NumberOfElems=size(ElemNodeConnectivity,1);
NumberOfNodes=size(X0all,1);
NumberOfAllDofs = DofsAtNode*NumberOfNodes; % Number of all dofs of structure

% Display the loaded matrix
%disp('Loaded Numeric Data from File:');
%disp(loadedData)

% Lets use function that creates DOFs connectivity
ElemDOFsConnectivity=ElemDOFsConnectivity_Q4(ElemNodeConnectivity,NumberOfElems);
%**

% Let's define bc vector (full of ones and zeros) to set up boundary
% conditions (now linear constraints only) easier
bc = logical(ones(1,NumberOfAllDofs));                  % 1 = no constraints, 0 = constraints
% TIP: Use FindGlobNodeID_v10 function to find node IDs for fixed
% coordinates and applied forces. Now they are not correctly dome
bcNodeIDsSymX = FindGlobNodeID_v10(X0all(:,1),0);                    % TIP: need to find nodes at symmetric line x=? Give value instead of 0
bcNodeIDsSymY = FindGlobNodeID_v10(X0all(:,2),0);                    % TIP: need to find nodes at symmetric line y=? Give value instead of 0
bcDofIDsSymX=DofID(DofsAtNode,bcNodeIDsSymX,[1:2]);             % TIP: Boundary condtions. Choose nodal displacements in x or y.  Now all fixed which is not correct! Remember symmetricity.
bcDofIDsSymY=DofID(DofsAtNode,bcNodeIDsSymY,[1:2]);             % TIP: Boundary condtions. Choose nodal displacements in x or y. Now all fixed which is not correct! Remember symmetricity. 
bc(bcDofIDsSymX)=0;
bc(bcDofIDsSymY)=0;
DOF=sum(bc);                                           % Number of unconstrained DOFs

AppliedForceNodeIDs=FindGlobNodeID_v10(X0all,100);                 % TIP: find node IDs where applied forces are active. Give coordinate x=? at the free end.  
AppliedForceDofIDs=DofID(DofsAtNode,AppliedForceNodeIDs,[1:2]);  % TIP: choose direction (DOFs where forces are applied). Now all [1,2] are chosen and it is not correct.


% Initialization
Kglob=zeros(NumberOfAllDofs,NumberOfAllDofs);

% reshape makes vector form from X0all that is needed to create X in assembling 
X0allvec=reshape(X0all',[],1);

   
% Assemble the global matrix. This works for all cases
for ii=1:NumberOfElems
    X=X0allvec(ElemDOFsConnectivity(ii,:));
    
    % full Gaussian int
    %Kloc=KlocQ4xi(X,E,nu,Lz);
    % Selective reduced int
    Kloc=KlocQ4psR(X,E,nu,Lz);              % Function of nodal coordinates (trasnformationm via isoparametric)    
    %ii
    %rankvalue(ii)=rank(Kloc); 
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
fglob(AppliedForceDofIDs)=F/length(AppliedForceDofIDs);  % distribute force at chosen dofs

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

PostProcessing_v10

