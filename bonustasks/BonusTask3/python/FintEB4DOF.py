import numpy as np

def f_intEB4DOF(L, v1, theta1, v2, theta2, E, Iz):
    F1 = (12*E*Iz*(v1 - v2))/L**3 + (6*E*Iz*(theta1 + theta2))/L**2
    F2 = (6*E*Iz*(v1 - v2))/L**2 + (2*E*Iz*(2*theta1 + theta2))/L
    F3 = (-12*E*Iz*(v1 - v2))/L**3 - (6*E*Iz*(theta1 + theta2))/L**2
    F4 = (6*E*Iz*(v1 - v2))/L**2 + (2*E*Iz*(theta1 + 2*theta2))/L
    return np.array([[F1], [F2], [F3], [F4]])
