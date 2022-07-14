# POTENTIAL FUNCTION:
import numpy as np
from matplotlib import pyplot as plt


def V_r(r, m):
    G = 6.6738480 * 10 ** (-11)  # Nm^2/Kg^2
    return -G * m / r


def V_no_m(r, _):
    return -1 / r


def V_eff(l,m,r,V):
    return l*l/(2*m*r*r)+V


def V_Yukawa(r, _):
    return -(np.exp(-r))/r


def plot_V_eff():
    r = np.linspace(0.01, 5, 100)
    r0 = [r[0], 0.65, 0.8, 1.25]  # r[0]
    m = 1
    v = 1
    l = r0*m*v
    for ii in l:
        y = V_eff(ii, m, r, V_no_m(r, m))
        plt.plot(r, y)
    plt.axis([0, 5, -1.5, 0.5])
    plt.axhline(y=0, color='black', linestyle='--')
    plt.legend(r0)
    plt.show()


if __name__=='__main__':
    plot_V_eff()
