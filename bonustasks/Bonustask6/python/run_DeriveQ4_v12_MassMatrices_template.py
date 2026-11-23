# BK10A6400 Basics of FE Analysis (FEMBasics2026)
# Teacher in charge: Marko Matikainen (MKM)
# Derivation of the bilinear plane element: mass matrices
# Code template for the Bonus Task 6 (Python version of run_DeriveQ4_v12_MassMatrices_template.m).
#
# Run this file first. It writes MlocQ4.py (consistent mass matrix) and
# MloclumpedQ4.py (lumped mass matrix) that are used by run_ModalAnalysis_Q4_v21.py.

import sympy as sp
from pythonFunction import pythonFunction

DOFs = 8         # element has 4 nodes 2 dofs per node (4x2=8)

xi, eta, Lx, Ly, Lz, E, nu = sp.symbols('xi eta Lx Ly Lz E nu')
u1, u2, u3, u4, v1, v2, v3, v4 = sp.symbols('u1 u2 u3 u4 v1 v2 v3 v4')
X1, X2, X3, X4, X5, X6, X7, X8 = sp.symbols('X1 X2 X3 X4 X5 X6 X7 X8')
rho = sp.symbols('rho')

# A vector of nodal displacements
uu = sp.Matrix([u1, v1, u2, v2, u3, v3, u4, v4])

# A vector of nodal position
X = sp.Matrix([X1, X2, X3, X4, X5, X6, X7, X8])

# A vector of nodal displacements separated into u and v
u = sp.Matrix([u1, u2, u3, u4])
v = sp.Matrix([v1, v2, v3, v4])

p = sp.Matrix([[1, xi, eta, xi*eta]])

AA = sp.Matrix([[1, -1, -1,  1],
                [1,  1, -1, -1],
                [1,  1,  1,  1],
                [1, -1,  1, -1]])

Nvec = p*AA.inv()

# Displacement field interpolation in x and y directions
uh = (Nvec*u)[0]
vh = (Nvec*v)[0]

# Interpolation for geometry (position)
Xinit = sp.Matrix([X1, X3, X5, X7])
Yinit = sp.Matrix([X2, X4, X6, X8])

Xh = (Nvec*Xinit)[0]
Yh = (Nvec*Yinit)[0]

Je = sp.zeros(2, 2)
Je[0, 0] = sp.diff(Xh, xi)
Je[0, 1] = sp.diff(Xh, eta)
Je[1, 0] = sp.diff(Yh, xi)
Je[1, 1] = sp.diff(Yh, eta)
detJe = sp.expand(Je.det())

# Consistent mass matrix (via shape functions)
# Shape function matrix in xi coordinate system (2 x 8)
N = sp.zeros(2, DOFs)
for ii in range(2):
    for jj in range(4):
        jj2 = ii + jj*2
        N[ii, jj2] = Nvec[jj]

NNdxideta = (N.T*N*detJe).applyfunc(sp.expand)
NNdeta = NNdxideta.applyfunc(lambda f: sp.integrate(f, (xi, -1, 1)))
NN_intdxideta = NNdeta.applyfunc(lambda f: sp.integrate(f, (eta, -1, 1)))
Mloc = rho*Lz*NN_intdxideta

pythonFunction(Mloc, 'MlocQ4', [X, rho, Lz])

# Lumped mass matrix
alpha = 1                                    # Tip: You need to find correct alpha...
Mlumpedtmp = alpha*rho*Lz*sp.eye(8)          # TIP: you need to think how much mass should
# be lumped at nodes...so this line must be modified by multiplying it by the
# factor alpha. Should be different than 1.
detJeint = sp.integrate(sp.integrate(detJe, (xi, -1, 1)), (eta, -1, 1))
Mloclumbed = Mlumpedtmp*detJeint
pythonFunction(Mloclumbed, 'MloclumpedQ4', [X, rho, Lz])

print('MlocQ4.py and MloclumpedQ4.py written.')
