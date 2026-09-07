import numpy as np

def kloc_rod(E, A, L):
    """Local 2x2 axial rod stiffness matrix."""
    K=(E * A / L) * np.array([[1.0, -1.0],
                            [-1.0,  1.0]])
    return K