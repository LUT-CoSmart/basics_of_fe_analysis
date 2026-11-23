# BK10A6400 Basics of FE Analysis (FEMBasics2026autumn)
# Teacher in charge: Marko Matikainen (MKM)
# Bonus task 6: Modal analysis of a beam-like plane structure with Q4 elements
# (Python version of run_ModalAnalysis_Q4_v21.m)
#
# Choose below:
#   IntType  : 1 = full Gaussian integration (KlocQ4xiFullGauss)
#              2 = selective reduced integration, SRI (KlocQ4psR)
#   MassType : 1 = consistent mass matrix (MlocQ4)
#              2 = lumped mass matrix (MloclumpedQ4, derive it first with
#                  run_DeriveQ4_v12_MassMatrices_template.py)
#   BC       : 1 = free structure (no supports)
#              2 = cantilever, all nodes at x = 0 clamped
#   ScaleMesh: mesh (4*ScaleMesh) x (1*ScaleMesh), ScaleMesh = 1/4 gives a 1x1 mesh
#
# Prepare this code so that you can quickly change the parameters, BC, MassType,
# IntType and the mesh: in the quiz, part A is computed in the class with new values.
#
# NOTE: node and DOF IDs are numbered from 1 as in the MATLAB code
# (RectangularMeshQ4, ElemDOFsConnectivity_Q4, DofID, FindGlobNodeID).
# They are converted to Python indices (-1) only when arrays are indexed.

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

from RectangularMeshQ4 import RectangularMeshQ4
from ElemDOFsConnectivity_Q4 import ElemDOFsConnectivity_Q4
from FindGlobNodeID import FindGlobNodeID
from DofID import DofID
from DrawMeshQ4 import DrawMeshQ4
from KlocQ4xiFullGauss import KlocQ4xiFullGauss
from KlocQ4psR import KlocQ4psR
from MlocQ4 import MlocQ4

#%% Choices
IntType   = 2       # 1 = full Gauss, 2 = SRI
MassType  = 1       # 1 = consistent, 2 = lumped
BC        = 1       # 1 = free, 2 = cantilever (clamped at x = 0)
ScaleMesh = 4       # 1 -> 4x1, 2 -> 8x2, 4 -> 16x4, 8 -> 32x8
nmodes    = 8       # number of computed modes
PlotModes = True    # draw the modes

if MassType == 2:
    from MloclumpedQ4 import MloclumpedQ4   # run run_DeriveQ4_v12_MassMatrices_template.py first

#%% Data (units m, N, kg)
E   = 1e9
nu  = 0.3
rho = 7850
L   = 2
H   = 0.4
W   = 0.4
Lz  = W

DofsAtNode = 2      # 2 dofs at every node
LocalSize  = 8      # element's size (4 nodes * 2 dofs)

#%% Mesh
if ScaleMesh < 1:
    nxElems = 1; nyElems = 1                   # 1x1 mesh
else:
    nxElems = int(4*ScaleMesh)                 # discretization in x-direction
    nyElems = int(1*ScaleMesh)                 # discretization in y-direction
NumberOfElems = nxElems*nyElems

X0all, ElemNodeConnectivity = RectangularMeshQ4(nxElems, nyElems, 0, L, 0, H)
NumberOfNodes   = X0all.shape[0]
NumberOfAllDofs = DofsAtNode*NumberOfNodes
ElemDOFsConnectivity = ElemDOFsConnectivity_Q4(ElemNodeConnectivity, NumberOfElems)

#%% Boundary conditions: bc vector (True = free DOF, False = constrained DOF)
bc = np.ones(NumberOfAllDofs, dtype=bool)
if BC == 2:
    ClampedNodeIDs, _ = FindGlobNodeID(X0all, 0)                  # nodes at x = 0
    bc[DofID(DofsAtNode, ClampedNodeIDs, [1, 2]) - 1] = False     # u = v = 0

#%% Assembly (sparse global matrices)
Kglob = sp.lil_matrix((NumberOfAllDofs, NumberOfAllDofs))
Mglob = sp.lil_matrix((NumberOfAllDofs, NumberOfAllDofs))
X0allvec = X0all.reshape(-1)               # [x1, y1, x2, y2, ...]

for ii in range(NumberOfElems):
    dofs = ElemDOFsConnectivity[ii, :] - 1   # Python indices
    X = X0allvec[dofs]

    if IntType == 1:
        Kloc = KlocQ4xiFullGauss(X, E, nu, Lz)    # full Gaussian integration
    else:
        Kloc = KlocQ4psR(X, E, nu, Lz)            # selective reduced integration

    if MassType == 1:
        Mloc = MlocQ4(X, rho, Lz)                 # consistent mass matrix
    else:
        Mloc = MloclumpedQ4(X, rho, Lz)           # lumped mass matrix

    Kglob[np.ix_(dofs, dofs)] += Kloc
    Mglob[np.ix_(dofs, dofs)] += Mloc

Kglob = Kglob.tocsc()
Mglob = Mglob.tocsc()
Kglobred = Kglob[bc, :][:, bc]
Mglobred = Mglob[bc, :][:, bc]

# Check: total mass in x-direction (sum of u-terms) should be rho*L*H*Lz
TotalMass = Mglob[0::2, :][:, 0::2].sum()
print(f'Total mass (x-direction): {TotalMass:.4f} kg, rho*L*H*Lz = {rho*L*H*Lz:.4f} kg')

#%% Eigenvalue problem (K - lambda*M)*phi = 0
# A free structure has a singular K (rigid body modes) => use a shift sigma = -1
nmodes = min(nmodes, int(bc.sum()))
if nmodes < bc.sum() - 1:                   # sparse solver for a few modes
    d, Eigenmodes = spla.eigsh(Kglobred, k=nmodes, M=Mglobred, sigma=-1.0, which='LM')
else:                                       # very small model: dense solver
    import scipy.linalg as sla
    d, Eigenmodes = sla.eigh(Kglobred.toarray(), Mglobred.toarray())
    d, Eigenmodes = d[:nmodes], Eigenmodes[:, :nmodes]
idx = np.argsort(d)
lam = d[idx]
Eigenmodes = Eigenmodes[:, idx]
Eigenfrequencies = np.sqrt(np.abs(lam))    # [rad/s], rigid body modes ~ 0

print(f'Mesh {nxElems}x{nyElems}, IntType {IntType}, MassType {MassType}, BC {BC}')
for k in range(nmodes):
    print(f'Mode {k+1:2d}: omega = {Eigenfrequencies[k]:10.3f} rad/s, f = {Eigenfrequencies[k]/(2*np.pi):9.3f} Hz')

# Modes with all DOFs (constrained DOFs = 0) for drawing
EigenmodesFull = np.zeros((NumberOfAllDofs, nmodes))
EigenmodesFull[bc, :] = Eigenmodes

#%% Post processing
if PlotModes:
    uglob = np.zeros(NumberOfAllDofs)
    DrawingLimits = [-0.5, 2.5, -0.5, 1.0]
    for imode in range(nmodes):
        plt.figure(imode + 1)
        # scale the mode so that the largest nodal displacement is 0.2 m
        ScaleDisp = 0.2/np.max(np.abs(EigenmodesFull[:, imode]))
        DrawMeshQ4(X0all, ElemNodeConnectivity, uglob, DofsAtNode, DrawingLimits, 'b')
        DrawMeshQ4(X0all, ElemNodeConnectivity, EigenmodesFull[:, imode]*ScaleDisp, DofsAtNode, DrawingLimits, 'r')
        plt.title(f'Mode {imode+1}: $\\omega$ = {Eigenfrequencies[imode]:.5g} rad/s', fontsize=14)
    plt.show()
