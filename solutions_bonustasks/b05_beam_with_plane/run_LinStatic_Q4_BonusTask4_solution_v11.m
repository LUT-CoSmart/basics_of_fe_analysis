% BK10A6400 Basics of FE Analysis (FEMBasics2023)
% Teacher in charge: Marko Matikainen (MKM)
% Template for the Bonus Task 4.
%
% Goal: The code solves global displacements at free end by using Q4 element. 
% Coded by MKM for student's usage in the FEMBasics2023 course


%          N4     N5    N6    
%          
%          \            |F
%          \o-----o-----o
%          \|     |     |
%          \|     |     |
%          \o-----o-----o
%          \
%         N1     N2    N3  


% So-called connectivity matrix
% Ele1 N1 N2 N5 N4  u1,v1,u2,v2,u5,v5,u4,v4
% Ele2 N2 N3 N6 N5  u2,v2,u3,v3,u6,v6,u5,v5,
% so nodes 2 and 5 are common nodes

% First matlab-based FE code for a simple beam structure analyzed by Q4
clear all;
close all;
clc;
%format long;
format shortG;

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
%nxElems=1;                                  % Discretization in x-direction
nxElems=2;
nyElems=1;                                  % Discretization in y-direction
NumberOfElems=nxElems*nyElems;
%NumberOfNodes=4;
NumberOfNodes=6;
NumberOfAllDofs = DofsAtNode*NumberOfNodes; % Number of all dofs of structure

Lx=L/nxElems;                               % Element's size in x-direction
Ly=H/nyElems;                               % Element's size in y-direction                                 
Lz=W;                                 

%ElemNodeConnectivity=[1,2,4,3];             % for one element. You need two lines for two elements...
% 
 ElemNodeConnectivity=[1,2,5,4;
                  2,3,6,5];



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

% % % One element
% % % % Load vector (in a global coordinate systems)
% fglob=zeros(8,1);
% fglob(4)=F/2;        % add loading
% fglob(8)=F/2;        % add loading
% 
% % Boundary conditions (in a global coordinate systems)
% Kglobred=zeros(4,4);
% fglobred=zeros(4,1);
% Kglobred=Kglob([3:4,7:8],[3:4,7:8]);
% fglobred=fglob([3:4,7:8]);
% 
% % Let's solve displacements uglob 
% uglobred=Kglobred\fglobred;
% 
% % % Let's gather all displacements (solved and fixed)
% ugloball=zeros(8,1);
% ugloball([3:4,7:8])=uglobred

% % Two elements
% % % % % Load vector (in a global coordinate systems)
fglob=zeros(12,1);
fglob(6)=F/2;
fglob(12)=F/2;

% Boundary conditions (in a global coordinate systems)
Kglobred=zeros(8,8);
fglobred=zeros(8,1);
Kglobred=Kglob([3:6,9:12],[3:6,9:12]);
fglobred=fglob([3:6,9:12]);

% Let's solve displacements uglob 
uglobred=Kglobred\fglobred


% % Let's gather all displacements (solved and fixed)
ugloball=zeros(12,1);
ugloball([3:6,9:12])=uglobred

detKglobred=det(Kglobred)
ugloball([5,6,11,12])

% Analytical solution
Iz=W*H^3/12;
A=W*H;
kappaxy=6/5;
nu=0.3;
G=E/(2*(1+nu));
deltamaxTimo=F*L^3/(3*E*Iz)+kappaxy*F*L/(G*A)



