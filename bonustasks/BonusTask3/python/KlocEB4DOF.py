import numpy as np

def klocEB4DOF(E, Iz, L):
    return np.array([
        [12*E*Iz/L**3, 6*E*Iz/L**2, -12*E*Iz/L**3, 6*E*Iz/L**2],
        [6*E*Iz/L**2, 4*E*Iz/L, -6*E*Iz/L**2, 2*E*Iz/L],
        [-12*E*Iz/L**3, -6*E*Iz/L**2, 12*E*Iz/L**3, -6*E*Iz/L**2],
        [6*E*Iz/L**2, 2*E*Iz/L, -6*E*Iz/L**2, 4*E*Iz/L]
    ])
