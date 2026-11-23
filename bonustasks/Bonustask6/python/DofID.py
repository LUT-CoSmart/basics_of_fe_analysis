# DofID.py
# Returns IDs of nodal coordinates (Python version of DofID.m, MKM).
#   nodes : list of chosen node IDs (from 1)
#   comps : list of chosen DOFs at node (1 = u, 2 = v)
# The returned DOF IDs start from 1 as in MATLAB (subtract 1 for Python indexing).

import numpy as np

def DofID(DofsAtNode, nodes, comps):
    ID = []
    for nn in np.atleast_1d(nodes):
        for c in np.atleast_1d(comps):
            ID.append(DofsAtNode*nn - DofsAtNode + c)
    return np.array(ID, dtype=int)
