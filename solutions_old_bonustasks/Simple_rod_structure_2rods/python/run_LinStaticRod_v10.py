# BK10A6400 Basics of FE Analysis
# Final solution code for Bonus Task 2
#
# Goal:
# Solve global displacements at node 2
# and member forces of a simple rod structure.
#
# Units: mm, N, radians


import numpy as np

from TRod import TRod
from KlocRod import KlocRod


# ==========================================================
# Rod 1
# ==========================================================

E1 = 210000
A1 = 100
L1 = 1000

alpha1 = -45 * np.pi / 180


# ==========================================================
# Rod 2
# ==========================================================

E2 = 210000
A2 = 50
L2 = 1000

alpha2 = -135 * np.pi / 180


# ==========================================================
# Force
# Force in vertical global direction at node 2
# ==========================================================

F = -100000


# ==========================================================
# Transformation into a global coordinate system
# ==========================================================

# Transformation matrix for rod 1
T1 = TRod(alpha1)

# Transformation matrix for rod 2
T2 = TRod(alpha2)


# ==========================================================
# Local elemental stiffness matrices
# ==========================================================

Kloc1 = KlocRod(E1, A1, L1)
Kloc2 = KlocRod(E2, A2, L2)


# ==========================================================
# Global elemental stiffness matrix for rod 1
# ==========================================================

Kglob1 = T1.T @ Kloc1 @ T1


# ==========================================================
# Global elemental stiffness matrix for rod 2
# ==========================================================

Kglob2 = T2.T @ Kloc2 @ T2


# ==========================================================
# Global stiffness matrix for whole structure
# 6 x 6
# ==========================================================

Kglob = np.zeros((6, 6))


# Substitute element 1 into global stiffness matrix
Kglob[0:4, 0:4] = Kglob1


# Substitute element 2 and add stiffnesses
Kglob[2:6, 2:6] += Kglob2


# ==========================================================
# Load vector in global coordinate system
# ==========================================================

fglob = np.zeros((6, 1))

# Force at DOF 4
fglob[3, 0] = F


# ==========================================================
# Boundary conditions
#
# u1, v1, u3, v3 are fixed
#
# Remaining DOFs:
# u2, v2
# ==========================================================

Kglobc = Kglob[2:4, 2:4]

fglobc = fglob[2:4]


# ==========================================================
# Check rank and determinant
# ==========================================================

print("Rank of Kglobc:")
print(np.linalg.matrix_rank(Kglobc))

print("\nDeterminant of Kglobc:")
print(np.linalg.det(Kglobc))


# ==========================================================
# Solve nodal displacement
# ==========================================================

uglobc = np.linalg.solve(Kglobc, fglobc)


# ==========================================================
# Collect all nodal displacements
# Fixed displacements are zero
# ==========================================================

ugloball = np.zeros((6, 1))

ugloball[2:4] = uglobc


# ==========================================================
# Elemental nodal displacement vectors
# ==========================================================

uglob1 = ugloball[0:4]

uglob2 = ugloball[2:6]


# ==========================================================
# Solve elemental axial displacements
# in local coordinate system
# ==========================================================

uloc1 = T1 @ uglob1

uloc2 = T2 @ uglob2


# ==========================================================
# Solve member forces
# f = K u
# ==========================================================

floc1 = Kloc1 @ uloc1

floc2 = Kloc2 @ uloc2


# ==========================================================
# Display results
# ==========================================================

print("\n========================================")
print("Determinant of Kglob")
print("========================================")

# Determinant of unconstrained global stiffness matrix
# Should be zero or very close to zero

print(np.linalg.det(Kglob))


print("\n========================================")
print("Determinant of Kglobc")
print("========================================")

# Determinant of constrained global stiffness matrix
# Should be > 0

print(np.linalg.det(Kglobc))


print("\n========================================")
print("Displacement at node 2")
print("========================================")

print(ugloball[2:4])


print("\n========================================")
print("Local elemental displacements")
print("========================================")

print("Element 1:")
print(uloc1)

print("\nElement 2:")
print(uloc2)


print("\n========================================")
print("Member forces")
print("========================================")

print("Element 1:")
print(floc1)

print("\nElement 2:")
print(floc2)