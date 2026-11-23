# DrawMeshQ4.py
# Draws the Q4 mesh from node to node (Python version of DrawMeshQ4.m).
#   X0   : nodal coordinates [x, y]
#   nloc : element nodal connectivity (node IDs from 1)
#   u    : global displacement vector (all DOFs)
#   color: matplotlib color, e.g. 'b' or 'r'

import numpy as np
import matplotlib.pyplot as plt
from DofID import DofID

def DrawMeshQ4(X0, nloc, u, DofsAtNode, DrawingLimits, color, ax=None):
    if ax is None:
        ax = plt.gca()
    NumberOfNodes = X0.shape[0]
    NumberOfElems = nloc.shape[0]

    # Makes an array [NumberOfNodes, 2] for drawing from the displacement vector
    loc = DofID(DofsAtNode, np.arange(1, NumberOfNodes + 1), [1, 2]) - 1   # Python indices
    Nu = np.reshape(u[loc], (NumberOfNodes, DofsAtNode))
    C = X0 + Nu        # updates X0 with nodal displacements

    # Instructions for drawing lines from node to node (local numbering)
    nlocdraw = [0, 1, 2, 3, 0]
    for ii in range(NumberOfElems):
        nodes = nloc[ii, nlocdraw] - 1                                   # Python indices
        ax.plot(C[nodes, 0], C[nodes, 1], color, linewidth=1.5)

    ax.set_aspect('equal')
    ax.set_xlim(DrawingLimits[0], DrawingLimits[1])
    ax.set_ylim(DrawingLimits[2], DrawingLimits[3])
    return ax
