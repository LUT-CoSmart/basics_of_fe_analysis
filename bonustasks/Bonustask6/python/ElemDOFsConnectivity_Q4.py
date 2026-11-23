# ElemDOFsConnectivity_Q4.py
# Makes a table of elements' nodal displacements (DOFs) connectivity from
# elements' node connectivity table (Python version of ElemDOFsConnectivity_Q4.m).
# DOF IDs start from 1 as in MATLAB: node k has DOFs 2k-1 (u) and 2k (v).

import numpy as np

def ElemDOFsConnectivity_Q4(ElemNodeConnectivity, NumberOfElems):
    DofsAtElem = 8
    ElemDofsConnectivity = np.zeros((NumberOfElems, DofsAtElem), dtype=int)
    for k in range(NumberOfElems):
        n1 = [ElemNodeConnectivity[k, 0]*2 - 1, ElemNodeConnectivity[k, 0]*2]
        n2 = [ElemNodeConnectivity[k, 1]*2 - 1, ElemNodeConnectivity[k, 1]*2]
        n3 = [ElemNodeConnectivity[k, 2]*2 - 1, ElemNodeConnectivity[k, 2]*2]
        n4 = [ElemNodeConnectivity[k, 3]*2 - 1, ElemNodeConnectivity[k, 3]*2]
        ElemDofsConnectivity[k, :] = n1 + n2 + n3 + n4
    return ElemDofsConnectivity
