from sympy import ImmutableDenseMatrix

def KlocRod(A, E, L):

    return ImmutableDenseMatrix([[A*E/L, -A*E/L], [-A*E/L, A*E/L]])

