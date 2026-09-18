# save_as_pyfunc_rod.py

import sympy as sp
from sympy.printing.pycode import pycode


# ==========================================================
# Let's initialize variables
# ==========================================================

x, L, u1, u2, A, E = sp.symbols('x L u1 u2 A E')


# ==========================================================
# A vector of nodal displacements
# ==========================================================

uu = sp.Matrix([u1, u2])


# ==========================================================
# Matrix for determining the shape functions
# ==========================================================

AA = sp.Matrix([
    [1, 0],
    [1, L]
])


# ==========================================================
# Polynomials
# ==========================================================

p = sp.Matrix([1, x])


# ==========================================================
# Shape functions
# ==========================================================

AA_inv = sp.simplify(AA.inv())

N = p.T * AA_inv
N = sp.simplify(N)


# ==========================================================
# Write Shapef2NodeRod.py
# ==========================================================

expr_src = pycode(N)

module_src = f"""from sympy import ImmutableDenseMatrix

def Shapef2NodeRod(x, L):

    return {expr_src}

"""

with open("Shapef2NodeRod.py", "w", encoding="utf-8") as f:
    f.write(module_src)

print("Wrote Shapef2NodeRod.py")


# ==========================================================
# Displacement field
# ==========================================================

uh = N * uu


# ==========================================================
# Axial strain
# ==========================================================

Epsxx = sp.diff(uh, x)


# ==========================================================
# Internal strain energy for a rod element
# ==========================================================

Wintdx = sp.Rational(1, 2) * E * A * Epsxx**2


# ==========================================================
# Integrate over element's length
# ==========================================================

Wint = sp.integrate(Wintdx, (x, 0, L))
Wint = sp.simplify(Wint)

print("\nWint =")
sp.pprint(Wint)


# ==========================================================
# Fint = d Wint / d u
#
# This is an intermediate quantity.
# It is NOT written to a separate file.
# ==========================================================

DOFs = 2

Fint = sp.Matrix([
    sp.diff(Wint, uu[kk])
    for kk in range(DOFs)
])

Fint = sp.simplify(Fint)


# ==========================================================
# Kloc = d Fint / d u
# ==========================================================

Kloc = sp.Matrix([
    [
        sp.diff(Fint[ii], uu[kk])
        for kk in range(DOFs)
    ]
    for ii in range(DOFs)
])

Kloc = sp.simplify(Kloc)


# ==========================================================
# Write WintRod.py
# ==========================================================

expr_src = pycode(Wint)

module_src = f"""from sympy import ImmutableDenseMatrix

def WintRod(u1, u2, A, E, L):

    return {expr_src}

"""

with open("WintRod.py", "w", encoding="utf-8") as f:
    f.write(module_src)

print("Wrote WintRod.py")


# ==========================================================
# Write KlocRod.py
# ==========================================================

expr_src = pycode(Kloc)

module_src = f"""from sympy import ImmutableDenseMatrix

def KlocRod(A, E, L):

    return {expr_src}

"""

with open("KlocRod.py", "w", encoding="utf-8") as f:
    f.write(module_src)

print("Wrote KlocRod.py")