# BK10A6400 Basics of FE Analysis (FEMBasics2026)
# Teacher in charge: Marko Matikainen (MKM)
# Solution for the Bonus Task 3.
# 
# Goal: The code solves global displacements at node 2 and member forces of a simple beam structure.
# Coded by MUS for student's usage in the FEMBasics2026 course

#          N1     N2       N3
#          L1,I1     L2,I2
#          \      |F
#          \o-----o--------o
#          \  (1)     (2) o^o
# 
#  v1,theta1   v2,theta2  v3,theta3

# So-called connectivity matrix
# Ele1 N1 N2   v1,theta1,v2,theta2
# Ele2 N2 N3   v2,theta2,v3,theta3
# so node 2 is common node (and then displacements v2,theta2) for the both beam elements

# First python-based FE code for a simple beam structure

import numpy as np
import matplotlib.pyplot as plt
import time

from KlocEB4DOF import klocEB4DOF
from Shapef_EB4DOF import shapef_EB4DOF
from TEB4DOFs import TEB4DOFs


# Units are mm and N and radians

# Beam 1
E1 = 210000
I1 = 4762
L1 = 1000
alpha1 = 0


# Beam 2
E2 = 210000
I2 = 9524
L2 = 2000
alpha2 = 0


# Force (in vertical (global) direction at node 2)
F = -100


# Elemental stiffness matrices
Kloc1 = klocEB4DOF(E1, I1, L1)
Kloc2 = klocEB4DOF(E2, I2, L2)

# TPS: add more stiffness matrices here when you want to get solution with
# finer meshes

print("Kloc1 =")
print(Kloc1)

print("\nKloc2 =")
print(Kloc2)


# Transformation into a global coordinate system
# Transformation matrix for a beam 1
T1 = TEB4DOFs(alpha1)

# Transformation matrix for a beam 2
T2 = TEB4DOFs(alpha2)

# TIPS
# Compute T:s for new elemenets...if you wnat


# Global elemental stiffness matrix for a rod 1
Kglob1 = T1.T @ Kloc1 @ T1

# Global elemental stiffness matrix for a rod 2
Kglob2 = T2.T @ Kloc2 @ T2

# TIPs: you can do this for finre meshes also if you want..


# As you can see, local stiffness matrices have zero rows and columns due
# to transformation to global coordinate system and paraller coordinate
# systems so no point to make transformation because zero rows and columns
# need to eliminated. Note that this is just special case. For 6 DOFs beam
# element   later, use transformation.

# Let's use Klocs directly:
# TIP: need t change for the finer meshes
Kglob = np.zeros((6, 6))

Kglob[0:4, 0:4] = Kloc1
Kglob[2:6, 2:6] = Kglob[2:6, 2:6] + Kloc2


# If this previous assembling way in python is somehow difficult,
# you can always do:
# Kglob[0:4, 0:4] = Kloc1
# Kglob[2:6, 2:6] = Kloc2
# but now stiffnesses of a common node is not computed properly
# So lets do it explicitely
# Kglob[2, 2] = Kloc1[2, 2] + Kloc2[0, 0]
# Kglob[2, 3] = Kloc1[2, 3] + Kloc2[0, 1]
# Kglob[3, 2] = Kloc1[3, 2] + Kloc2[1, 0]
# Kglob[3, 3] = Kloc1[3, 3] + Kloc2[1, 1]

print("\nKglob =")
print(Kglob)


# TIP: need t change for the finer meshes
# Load vector (in a global coordinate systems)
fglob = np.zeros(6)

fglob[2] = F


# TIP: need t change for the finer meshes
# Boundary conditions (in a global coordinate systems)
# v1=0,theta1=0,v3=0
# Remaining (free) dofs (nodal displacements) are related to indeces 2,3,5 (0,1,4 (v1,theta1,v3) are fixed)
# Remanininf DOFs of system = 6 -3 (number of constraints) = 3

Kglobred = np.zeros((3, 3))
fglobred = np.zeros(3)

free_dofs = [2, 3, 5]

Kglobred = Kglob[np.ix_(free_dofs, free_dofs)]
fglobred = fglob[free_dofs]


# Let's solve determinant
det_Kglobred = np.linalg.det(Kglobred)

print("\ndet(Kglobred) =")
print(det_Kglobred)


# Let's solve displacements uglob
uglobred = np.linalg.solve(Kglobred, fglobred)

print("\nuglobred =")
print(uglobred)


# Next lines are just for see why it is beneficial to use
# backslash for solvig a system of linear equations uglobred=Kglobred^-1*fglobred
# Comment lines out and see computations time for these different opearations

# tic = time.perf_counter()
# uglobred = np.linalg.inv(Kglobred) @ fglobred
# toc = time.perf_counter()
# print("\ntime for different methods =")
# print(toc - tic)

# tic = time.perf_counter()
# uglobred = np.linalg.matrix_power(Kglobred, -1) @ fglobred
# toc = time.perf_counter()
# print(toc - tic)

# tic = time.perf_counter()
# uglobred = np.linalg.solve(Kglobred, fglobred)
# toc = time.perf_counter()
# print(toc - tic)


# TIP: need t change for the finer meshes
# Let's gather all displacements (solved and fixed)
ugloball = np.zeros(6)

ugloball[free_dofs] = uglobred

print("\nugloball =")
print(ugloball)


# TIP: need t change for the finer meshes
# Let's gather elemental displacement vectors
ug1 = ugloball[0:4]
ug2 = ugloball[2:6]


# Now no need to make transformation
# uloc1=T1@ug1;
# uloc2=T2@ug2;

uloc1 = ug1
uloc2 = ug2

print("\nuloc1 =")
print(uloc1)

print("\nuloc2 =")
print(uloc2)


# TIPS
# uloc3=ug3
# ...


# Let's solve member forces (in a local coordinate system)
Floc1 = Kloc1 @ uloc1
Floc2 = Kloc2 @ uloc2

print("\nFloc1 =")
print(Floc1)

print("\nFloc2 =")
print(Floc2)


# TIPS
# Floc3=Kloc3 @ uloc3
# ...


# Post processing
# ###############################################################
# ###############################################################

# Simple visualization using only nodal displacements

scaleFactor = 20   # visualization scale factor


# Nodal x-coordinates
xnod = np.array([0, L1, L1 + L2])


# Nodal vertical displacements
vnod = np.array([
    ugloball[0],
    ugloball[2],
    ugloball[4]
])


plt.figure(1)
plt.hold = True


# Undeformed beam
plt.plot(
    xnod,
    [0, 0, 0],
    'k--',
    linewidth=1.5,
    label='Undeformed beam'
)


# Deformed beam using only nodal displacements
plt.plot(
    xnod,
    scaleFactor * vnod,
    '-o',
    linewidth=2,
    label='Deformed beam'
)


plt.grid(True)
plt.box(True)
plt.axis('equal')

plt.xlabel('x [mm]')
plt.ylabel('Scaled vertical displacement')

plt.title('Finite element solution')

plt.legend(
    ['Undeformed beam', 'Deformed beam'],
    loc='best'
)

print(f"Visualization scale factor = {scaleFactor:g}")


# Visualization using shape functions

scaleFactor = 50
nplot = 50


# ----------------------------------
# Element 1
# ----------------------------------

x1loc = np.linspace(0, L1, nplot)
vplot1 = np.zeros(len(x1loc))

for ii in range(len(x1loc)):
    N = shapef_EB4DOF(x1loc[ii], L1)
    vplot1[ii] = np.array(N, dtype=float) @ uloc1

x1glob = x1loc


# ----------------------------------
# Element 2
# ----------------------------------

x2loc = np.linspace(0, L2, nplot)
vplot2 = np.zeros(len(x2loc))

for ii in range(len(x2loc)):
    N = shapef_EB4DOF(x2loc[ii], L2)
    vplot2[ii] = np.array(N, dtype=float) @ uloc2

x2glob = L1 + x2loc


# ----------------------------------
# Plot
# ----------------------------------

plt.figure(2)


# Undeformed beam
plt.plot(
    [0, L1, L1 + L2],
    [0, 0, 0],
    'k--',
    label='Undeformed'
)


# FE displacement field
plt.plot(
    x1glob,
    scaleFactor * vplot1,
    linewidth=2,
    label='Element 1'
)

plt.plot(
    x2glob,
    scaleFactor * vplot2,
    linewidth=2,
    label='Element 2'
)


plt.grid(True)
plt.axis('equal')

plt.xlabel('x [mm]')
plt.ylabel('Scaled displacement [mm]')

plt.legend(
    ['Undeformed', 'Element 1', 'Element 2']
)


# Show figures
plt.show()