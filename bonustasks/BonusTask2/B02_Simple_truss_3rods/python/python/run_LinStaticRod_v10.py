
# % BK10A6400 Basics of FE Analysis (FEMBasics2026autumn)
# % Teacher in charge: Marko Matikainen (MKM)
# %
# % Final solution code for the Bonus Task 2.
# %
# % Goal: The code solves global displacements at node 2
# % and member forces of a simple rod structure.
#
# Units: mm, N, radians

import numpy as np
from TRod import TRod
from KlocRod import KlocRod

np.set_printoptions(precision=6, suppress=True)


# ============================================================
# ROD 1
# ============================================================

E1 = 210000
A1 = 100
L1 = 1000
alpha1 = -45 * np.pi / 180

# ============================================================
# ROD 2
# ============================================================

E2 = 210000
A2 = 50
L2 = 1000
alpha2 = -135 * np.pi / 180

# ============================================================
# ROD 3
# ============================================================

E3 = 210000
A3 = 10
L3 = np.cos(np.pi / 4) * 1000
alpha3 = 0

# ============================================================
# FORCE
# ============================================================

# Vertical global force at node 2
F = -100000


# ============================================================
# TRANSFORMATION MATRICES
# ============================================================

T1 = TRod(alpha1)
T2 = TRod(alpha2)
T3 = TRod(alpha3)


# ============================================================
# LOCAL ELEMENTAL STIFFNESS MATRICES
# ============================================================

Kloc1 = KlocRod(E1, A1, L1)
Kloc2 = KlocRod(E2, A2, L2)
Kloc3 = KlocRod(E3, A3, L3)


# ============================================================
# GLOBAL ELEMENTAL STIFFNESS MATRICES
# ============================================================

Kglob1 = T1.T @ Kloc1 @ T1
Kglob2 = T2.T @ Kloc2 @ T2

# Rod 3 is horizontal, so its transformation is effectively
# already included in the local stiffness matrix.
Kglob3 = Kloc3


# ============================================================
# GLOBAL STIFFNESS MATRIX
# ============================================================

# 3 nodes × 2 DOF/node = 6 DOFs
Kglob = np.zeros((6, 6))

# Element 1: nodes 1-2
Kglob[0:4, 0:4] = Kglob1

# Element 2: nodes 2-3
Kglob[2:6, 2:6] += Kglob2

# NOTE:
# Your original MATLAB code calculates Kglob3 but does NOT
# actually assemble it into Kglob.
#
# To include Rod 3, its contribution must be assembled here.
# Assuming Rod 3 connects node 2 to node 3:
#
# Kglob[2:6, 2:6] += Kglob3


# ============================================================
# LOAD VECTOR
# ============================================================

fglob = np.zeros(6)

# Vertical force at node 2
# DOF ordering:
# [u1, v1, u2, v2, u3, v3]
fglob[3] = F


# ============================================================
# BOUNDARY CONDITIONS
# ============================================================

# u1, v1, u3, v3 are fixed.
#
# Free DOFs:
# u2, v2
#

Kglobc = Kglob[2:4, 2:4]
fglobc = fglob[2:4]


# ============================================================
# CHECK MATRIX RANK AND DETERMINANT
# ============================================================

print("Rank of Kglobc:")
print(np.linalg.matrix_rank(Kglobc))

print("\nDeterminant of Kglobc:")
print(np.linalg.det(Kglobc))


# ============================================================
# SOLVE NODAL DISPLACEMENTS
# ============================================================

uglobc = np.linalg.solve(Kglobc, fglobc)


# ============================================================
# COLLECT ALL NODAL DISPLACEMENTS
# ============================================================

ugloball = np.zeros(6)

# Displacements of node 2
ugloball[2:4] = uglobc


# ============================================================
# ELEMENTAL GLOBAL DISPLACEMENT VECTORS
# ============================================================

# Element 1: nodes 1-2
uglob1 = ugloball[0:4]

# Element 2: nodes 2-3
uglob2 = ugloball[2:6]


# ============================================================
# LOCAL ELEMENTAL DISPLACEMENTS
# ============================================================

uloc1 = T1 @ uglob1
uloc2 = T2 @ uglob2


# ============================================================
# MEMBER FORCES
# ============================================================

floc1 = Kloc1 @ uloc1
floc2 = Kloc2 @ uloc2


# ============================================================
# OUTPUT
# ============================================================

print("\nDeterminant of Kglob:")
print(np.linalg.det(Kglob))

print("\nDeterminant of Kglobc:")
print(np.linalg.det(Kglobc))

print("\nDisplacement at node 2:")
print(ugloball[2:4])

print("\nLocal elemental displacements - Element 1:")
print(uloc1)

print("\nLocal elemental displacements - Element 2:")
print(uloc2)

print("\nMember forces - Element 1:")
print(floc1)

print("\nMember forces - Element 2:")
print(floc2)