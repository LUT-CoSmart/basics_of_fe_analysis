# BK10A6400 Basics of FE Analysis (FEMBasics2025)
# Teacher in charge: Marko Matikainen
# Python/SymPy version of the MATLAB bonus task template
#
# Derives stiffness matrix for a two-node linear rod (bar/truss) element.

from sympy import symbols, Matrix, diff, integrate, simplify, lambdify, pprint #if it doesnot work run : pip install sympy
import numpy as np

# Symbols
x, L, A, E = symbols('x L A E', positive=True)
u1, u2 = symbols('u1 u2')
uu = Matrix([u1, u2])

# Polynomial basis p = [1, x]^T
p = Matrix([1, x])

# Build AA by evaluating p at the two nodes: x=0 and x=L
# TIP: Empty A, add components!
AA = Matrix([])

# Shape functions N = p^T * AA^{-1}
N = (p.T * AA.inv())              # row vector [N1, N2]
N = simplify(N)

# Displacement interpolation: u_h(x) = N * uu
uh = simplify(N * uu)             # scalar

# Axial strain: ε_xx = du_h/dx
Epsxx = simplify(diff(uh, x))     # should become (u2 - u1)/L

# Internal strain energy density per unit length:
# W_int,dx = 1/2 * A * E * ε_xx^2
Wintdx = simplify(0.5 * A * E * Epsxx**2)

# Internal strain energy over element (x in [0, L])
Wint = simplify(integrate(Wintdx, (x, 0, L)))

DOFs = 2

# Fint(kk) = dWint/duu(kk)
Fint = Matrix.zeros(DOFs, 1)
for kk in range(DOFs):
    Fint[kk, 0] = diff(Wint, uu[kk])

# Kloc(ii,kk) = dFint(ii)/duu(kk)
Kloc = Matrix.zeros(DOFs, DOFs)
for ii in range(DOFs):
    for kk in range(DOFs):
        Kloc[ii, kk] = diff(Fint[ii], uu[kk])

Kloc = simplify(Kloc)

# === Pretty printing of the symbolic results (optional) ===


print("Shape functions N = [N1, N2]:")
pprint(N)
print("\nInterpolated displacement u_h(x):")
pprint(uh)
print("\nAxial strain Epsxx:")
pprint(Epsxx)
print("\nInternal energy density Wintdx:")
pprint(Wintdx)
print("\nInternal energy Wint:")
pprint(Wint)
print("\nLocal stiffness matrix Kloc:")
pprint(Kloc)

# === Numeric callables (optional helpers) ===
#   shapef_2node_rod(x_val, L_val) -> [N1, N2]
shapef_2node_rod = lambdify((x, L), N, 'numpy')

#ExampleShapeFunction=shapef_2node_rod(11,2)

#   kloc_rod(A_val, E_val, L_val) -> 2x2 matrix
kloc_rod = lambdify((A, E, L), Kloc, 'numpy')


