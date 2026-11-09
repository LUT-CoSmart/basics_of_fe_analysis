import numpy as np

def TEB4DOFs(alpha):
    """
    Transformation matrix for a 2D beam element.

    For this particular problem alpha = 0, so the
    transformation matrix is simply the identity matrix.
    """
    c = np.cos(alpha)
    s = np.sin(alpha)

    # For a beam with vertical displacement and rotation
    # DOFs, the transformation is represented as:
    T = np.array([
        [ c,  0,  0,  0],
        [ 0,  1,  0,  0],
        [ 0,  0,  c,  0],
        [ 0,  0,  0,  1]
    ])

    return T