# RectangularMeshQ4.py
# Makes a rectangular mesh with four-node plane elements Q4
# (Python version of RectangularMeshQ4.m, MKM 21032022: modified for the linear FE course)
#
#   X0all, NodeConnect = RectangularMeshQ4(n, m, a1, a2, b1, b2)
#   n, m   : number of elements in x- and y-direction
#   a1, a2 : x-limits, b1, b2 : y-limits
#   X0all  : nodal coordinates, one row per node [x, y]
#   NodeConnect : element nodal connectivity (node IDs start from 1 as in MATLAB)

import numpy as np
from geospace import geospace

def RectangularMeshQ4(n, m, a1, a2, b1, b2):
    # Geospace for spacing
    xk = geospace(a1, a2, n + 1, 1)
    yk = geospace(b1, b2, m + 1, 1)

    X0all = []
    for k in range(m + 1):
        X0allk = np.column_stack([xk, yk[k]*np.ones(n + 1)])
        X0all.append(X0allk)
    X0all = np.vstack(X0all)

    # generate element connectivity (node IDs from 1)
    NodeConnect = []
    for l in range(1, m + 1):
        for k in range(1, n + 1):
            NodeConnectk = np.array([k, k + 1, k + 2 + n, n + k + 1]) + (l - 1)*(n + 1)
            NodeConnect.append(NodeConnectk)
    NodeConnect = np.array(NodeConnect, dtype=int)
    return X0all, NodeConnect
