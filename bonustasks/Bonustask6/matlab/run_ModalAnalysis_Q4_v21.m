% BK10A6400 Basics of FE Analysis (FEMBasics2026autumn)
% Teacher in charge: Marko Matikainen (MKM)
% Bonus task 6: Modal analysis of a beam-like plane structure with Q4 elements
%
% Choose below:
%   IntType  : 1 = full Gaussian integration (KlocQ4xiFullGauss)
%              2 = selective reduced integration, SRI (KlocQ4psR)
%   MassType : 1 = consistent mass matrix (MlocQ4)
%              2 = lumped mass matrix (MloclumpedQ4, derive it first with
%                  run_DeriveQ4_v12_MassMatrices_template.m)
%   BC       : 1 = free structure (no supports)
%              2 = cantilever, all nodes at x = 0 clamped
%   ScaleMesh: mesh (4*ScaleMesh) x (1*ScaleMesh), ScaleMesh = 1/4 gives a 1x1 mesh
%
% Prepare this code so that you can quickly change H, BC, MassType, IntType
% and the mesh: in the quiz, part A is computed in the class with new values.

clear all;
close all;
clc;
format long;

%% Choices
IntType   = 2;      % 1 = full Gauss, 2 = SRI
MassType  = 1;      % 1 = consistent, 2 = lumped
BC        = 1;      % 1 = free, 2 = cantilever (clamped at x = 0)
ScaleMesh = 1;      % 1 -> 4x1, 2 -> 8x2, 4 -> 16x4, 8 -> 32x8
nmodes    = 8;      % number of computed modes
PlotModes = true;   % draw the modes

%% Data (units m, N, kg)
E   = 1e9;
nu  = 0.3;
rho = 7850;
L   = 2;
H   = 0.4;
W   = 0.4;
Lz  = W;

DofsAtNode = 2;      % 2 dofs at every node
LocalSize  = 8;      % element's size (4 nodes * 2 dofs)

%% Mesh
if ScaleMesh < 1
    nxElems = 1;  nyElems = 1;                 % 1x1 mesh
else
    nxElems = 4*ScaleMesh;                     % discretization in x-direction
    nyElems = 1*ScaleMesh;                     % discretization in y-direction
end
NumberOfElems = nxElems*nyElems;

[X0all,ElemNodeConnectivity] = RectangularMeshQ4(nxElems,nyElems,0,L,0,H);
NumberOfNodes   = size(X0all,1);
NumberOfAllDofs = DofsAtNode*NumberOfNodes;
ElemDOFsConnectivity = ElemDOFsConnectivity_Q4(ElemNodeConnectivity,NumberOfElems);

%% Boundary conditions: bc vector (1 = free DOF, 0 = constrained DOF)
bc = logical(ones(1,NumberOfAllDofs));
if BC == 2
    ClampedNodeIDs = FindGlobNodeID(X0all,0);              % nodes at x = 0
    bc(DofID(DofsAtNode,ClampedNodeIDs,[1,2])) = 0;        % u = v = 0
end

%% Assembly (sparse global matrices)
Kglob = sparse(NumberOfAllDofs,NumberOfAllDofs);
Mglob = sparse(NumberOfAllDofs,NumberOfAllDofs);
X0allvec = reshape(X0all',[],1);

for ii = 1:NumberOfElems
    dofs = ElemDOFsConnectivity(ii,:);
    X    = X0allvec(dofs);

    if IntType == 1
        Kloc = KlocQ4xiFullGauss(X,E,nu,Lz);       % full Gaussian integration
    else
        Kloc = KlocQ4psR(X,E,nu,Lz);               % selective reduced integration
    end

    if MassType == 1
        Mloc = MlocQ4(X,rho,Lz);                   % consistent mass matrix
    else
        Mloc = MloclumpedQ4(X,rho,Lz);             % lumped mass matrix
    end

    Kglob(dofs,dofs) = Kglob(dofs,dofs) + Kloc;
    Mglob(dofs,dofs) = Mglob(dofs,dofs) + Mloc;
end

Kglobred = Kglob(bc,bc);
Mglobred = Mglob(bc,bc);

% Check: total mass in x-direction (sum of u-terms) should be rho*L*H*Lz
TotalMass = full(sum(sum(Mglob(1:2:end,1:2:end))));
fprintf('Total mass (x-direction): %.4f kg, rho*L*H*Lz = %.4f kg\n', TotalMass, rho*L*H*Lz);

%% Eigenvalue problem (K - lambda*M)*phi = 0
% A free structure has a singular K (rigid body modes) => use a shift sigma = -1
nmodes = min(nmodes, sum(bc));
if nmodes < sum(bc)-1
    [Eigenmodes,d] = eigs(Kglobred,Mglobred,nmodes,-1);    % sparse solver, a few modes
else
    [Eigenmodes,d] = eig(full(Kglobred),full(Mglobred));   % very small model (e.g. 1x1 mesh)
    Eigenmodes = Eigenmodes(:,1:nmodes);  d = d(1:nmodes,1:nmodes);
end
[lambda,idx]   = sort(real(diag(d)));
Eigenmodes     = Eigenmodes(:,idx);
Eigenfrequencies = sqrt(abs(lambda));      % [rad/s], rigid body modes ~ 0

fprintf('Mesh %dx%d, IntType %d, MassType %d, BC %d\n', nxElems, nyElems, IntType, MassType, BC);
for k = 1:nmodes
    fprintf('Mode %2d: omega = %10.3f rad/s, f = %9.3f Hz\n', k, Eigenfrequencies(k), Eigenfrequencies(k)/(2*pi));
end

% Modes with all DOFs (constrained DOFs = 0) for drawing
EigenmodesFull = zeros(NumberOfAllDofs,nmodes);
EigenmodesFull(bc,:) = Eigenmodes;

%% Post processing
if PlotModes
    uglob = zeros(NumberOfAllDofs,1);
    DrawingLimits = [-0.5 2.5 -0.5 1.0 -0.1 0.1];
    for imode = 1:nmodes
        figure(imode)
        % scale the mode so that the largest nodal displacement is 0.2 m
        ScaleDisp = 0.2/max(abs(EigenmodesFull(:,imode)));
        DrawMeshQ4(X0all,ElemNodeConnectivity,uglob,DofsAtNode,DrawingLimits,'b')
        DrawMeshQ4(X0all,ElemNodeConnectivity,EigenmodesFull(:,imode)*ScaleDisp,DofsAtNode,DrawingLimits,'r')
        title(['Mode ' int2str(imode) ': \omega = ' num2str(Eigenfrequencies(imode),5) ' rad/s'], ...
              'FontName','Times New Roman','FontSize',16);
    end
end
