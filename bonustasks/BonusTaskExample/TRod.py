import numpy as np

def t_rod(alpha):
    """
    Transformation matrix for a 2-node rod element.
    
    Parameters
    ----------
    alpha : float
        Angle in radians (measured from +x axis, CCW positive).
    
    Returns
    -------
    T : (2,4) ndarray
        Transformation matrix.
    """
    c = np.cos(alpha)
    s = np.sin(alpha)
    T = np.array([
        [c, s, 0.0, 0.0],
        [0.0, 0.0, c,  s ]
    ])
    return T
