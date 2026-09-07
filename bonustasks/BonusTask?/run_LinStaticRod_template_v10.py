import numpy as np
from KlocRod import kloc_rod
from TRod import t_rod


# Node and elemental numbering (ID) for a structure:
#           o N1
#            \
#       Elem1 \   
#              o  N2  
#       Elem2 /|
#            /  F 
#           o N3

# Units: mm and N; angles in radians
# Rod 1
E1 = 210000.0
A1 = 100.0
L1 = 1000.0
alpha1 = ??

# Rod 2
E2 = 210000.0
A2 = 50.0
L2 = 1000.0
alpha2 = ??

# Force: vertical at node 2 (negative y is downward)
F = -100000.0

# Transformations
T1 = t_rod(alpha1)
T2 = t_rod(alpha2)

# Local stiffness
Kloc1 = kloc_rod(E1, A1, L1)
Kloc2 = kloc_rod(E2, A2, L2)

# Global element stiffness
Kglob1 = T1.T @ Kloc1 @ T1
Kglob2 = T2.T @ Kloc2 @ T2

# Assemble global 6×6 (DOF: [u1, v1, u2, v2, u3, v3])
Kglob = np.zeros((6, 6), dtype=float)
Kglob[0:4, 0:4] = Kglob1          # elem 1: nodes 1–2
Kglob[2:6, 2:6] = Kglob[2:6, 2:6]+ Kglob2          # elem 2: nodes 2–3

# Global load vector
fglob = np.zeros((6,1),dtype=float)
fglob[?] = ??  # at v2

# Apply BCs: u1=v1=u3=v3=0 -> free DOFs: [u2, v2] -> indices [2,4]
Kglobc = Kglob[?,?]
fglobc = fglob[?:?]

# Checks
det_full = np.linalg.det(Kglob)
det_red = np.linalg.det(Kglobc)
rank_red = np.linalg.matrix_rank(Kglobc)

# Solve reduced system
uglobc = np.linalg.solve(Kglobc, fglobc)

# Collect all displacements
ugloball = np.zeros((6,1),dtype=float)
ugloball[?:?] = uglobc

# Elemental global displacement vectors
uglob1 = ugloball[?:?]   # [u1, v1, u2, v2]
uglob2 = ugloball[?:?]   # [u2, v2, u3, v3]

# Transform to local
uloc1 = ??
uloc2 = ??

# Local member forces f = Kloc * u_loc
floc1 = ??
floc2 = ??

# Output (mirrors MATLAB prints)
print('Determinant of Kglob')
print(det_full:.6g)
print('Determinant of Kglobc')
print(det_red:.6g)
print('Displacement at node 2')
print(ugloball[2:4])
print('Local elemental displacements')
print('uloc1 =', uloc1)
print('uloc2 =', uloc2)
print('Member forces (N)')
print('floc1 =', floc1)
print('floc2 =', floc2)