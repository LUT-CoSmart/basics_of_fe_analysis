# FindGlobNodeID.py
# Returns global nodal IDs (from 1) with a given coordinate in x-direction and
# the number of IDs (Python version of FindGlobNodeID.m, MKM).

import numpy as np

def FindGlobNodeID(X0, x):
    tol = 2*np.sqrt(np.finfo(float).eps)
    NodeID = [ii + 1 for ii in range(X0.shape[0]) if abs(X0[ii, 0] - x) < tol]
    NumberOfIDs = len(NodeID)
    return np.array(NodeID, dtype=int), NumberOfIDs
