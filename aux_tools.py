import numpy as np


def polar2cartesian(r, phi, psi):
    x = r * np.cos(phi) * np.cos(psi)
    y = r * np.sin(phi) * np.cos(psi)
    z = r * np.sin(psi)

    return [x, y, z]