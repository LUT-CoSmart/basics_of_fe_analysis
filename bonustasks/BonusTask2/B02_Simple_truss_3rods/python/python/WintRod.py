from sympy import ImmutableDenseMatrix

def WintRod(u1, u2, A, E, L):

    return ImmutableDenseMatrix([[(1/2)*A*E*(u1 - u2)**2/L]])

