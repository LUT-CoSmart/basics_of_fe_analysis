from sympy import ImmutableDenseMatrix

def Shapef2NodeRod(x, L):

    return ImmutableDenseMatrix([[(L - x)/L, x/L]])

