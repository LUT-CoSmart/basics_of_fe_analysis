# gauleg2.py
# Gives the points and weights for Gaussian quadrature (Python version of gauleg2.m).
#
#   xx, w = gauleg2(aa, bb, n) returns the points in area [aa, bb] and the
#   weights (in area [-1, 1]) for n-point Gaussian quadrature.

import numpy as np

def gauleg2(aa, bb, n):
    x = np.zeros(n); w = np.zeros(n)
    m = (n + 1)//2
    xm = 0.0; xl = 1.0
    for i in range(1, m + 1):
        z = np.cos(np.pi*(i - 0.25)/(n + 0.5))
        while True:
            p1 = 1.0; p2 = 0.0
            for j in range(1, n + 1):
                p3 = p2; p2 = p1; p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j
            pp = n*(z*p1 - p2)/(z*z - 1.0)
            z1 = z; z = z1 - p1/pp
            if abs(z - z1) < np.finfo(float).eps:
                break
        x[i-1] = xm - xl*z; x[n-i] = xm + xl*z
        w[i-1] = 2.0*xl/((1.0 - z*z)*pp*pp); w[n-i] = w[i-1]
    xx = (bb - aa)/2*x + (aa + bb)/2
    return xx, w
