# geospace.py
# Geometrically spaced vector (Python version of geospace.m).
#
# geospace(X1, X2, N, G) generates N points between X1 and X2 with a spacing
# defined by G. Each step is G times bigger than the previous step.
# G = 1 gives equally spaced points (as linspace).

import numpy as np

def geospace(X1, X2, N=100, G=1):
    if N < 2:                      # too few points specified
        return np.array([X2], dtype=float)
    if G == 1:
        return np.linspace(X1, X2, N)
    # Calculate geometrically spaced points
    return X1 + (np.cumprod(G*np.ones(N)) - G)/(G**N - G)*(X2 - X1)
