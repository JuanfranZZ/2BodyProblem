# POTENTIAL FUNCTION:
def V_real(r, m):
    G = 6.6738480 * 10 ** (-11)  # Nm^2/Kg^2
    return -G * m / r


def V_adim(r, m):
    return -1 / r