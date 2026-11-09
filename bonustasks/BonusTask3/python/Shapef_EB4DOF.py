import numpy as np

def shapef_EB4DOF(x, L):
    xi = x / L
    N1 = 1 - 3*xi**2 + 2*xi**3
    N2 = L*(xi - 2*xi**2 + xi**3)
    N3 = 3*xi**2 - 2*xi**3  
    N4 = L*(-xi**2 + xi**3)
    return np.array([N1, N2, N3, N4])
