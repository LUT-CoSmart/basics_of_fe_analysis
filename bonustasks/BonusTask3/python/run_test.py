# % BK10A6400 Basics of FE Analysis (FEMBasics2023autumn)
# % Teacher in charge: Marko Matikainen (MKM)
# % Test for the Bonus Task 3.
#
# % Let's test developed functions for the internal forces.

import numpy as np
import matplotlib.pyplot as plt

from Shapef_EB4DOF import shapef_EB4DOF
from FintEB4DOF import f_intEB4DOF


# clear all;
# close all;
# clc;


# % Initial values
E = 210000
Iz = 4762
L = 1000


# % Give nodal displacements v1, theta1, v2, theta2:
uu = np.array([0, 0, -L/1000, 0], dtype=float)


# % Let's visualize displacement field
x = np.arange(0, L + 100, 100)

uh = np.zeros(len(x))

for ii in range(len(x)):
    N = shapef_EB4DOF(x[ii], L)
    uh[ii] = np.array(N, dtype=float) @ uu


# % Approximate displacement uh at x={500,750} with given u
uh500 = np.array(shapef_EB4DOF(500, L), dtype=float) @ uu
uh750 = np.array(shapef_EB4DOF(750, L), dtype=float) @ uu

print("uh500 =", uh500)
print("uh750 =", uh750)


# % Plots dislacement field
plt.figure(1)
plt.plot(x, uh, 'r-')
plt.plot([750, 500], [uh750, uh500], 'b*')
plt.xlabel('Longitudinal coordinate x [mm]')
plt.ylabel('Displacment [mm]')


# % Always check your internal forces that it does not produce any forces
# % with zero displacements.
Fint = np.array(
    f_intEB4DOF(L, uu[0], uu[1], uu[2], uu[3], E, Iz),
    dtype=float
)

plt.show()