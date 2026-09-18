import numpy as np


def TRod(alpha):
    """
    Transformation matrix for a two-node rod element.
    """

    T = np.array([
        [np.cos(alpha), np.sin(alpha), 0, 0],
        [0, 0, np.cos(alpha), np.sin(alpha)]
    ])

    return T