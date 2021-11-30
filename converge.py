import numpy as np


def check(o_omega,n_omega,criteria):
    convergence = np.subtract(o_omega,n_omega)
    if (np.amax(convergence,axis = None) <= criteria):
        return True
    else:
        return False