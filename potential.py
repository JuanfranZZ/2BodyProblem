# POTENTIAL FUNCTION:
import numpy as np
from matplotlib import pyplot as plt


def V_r(r, m):

    return -m / r


def V_no_m(r, _):
    return -1 / r


def V_eff(l, m, r, V):
    return l*l/(2*m*r*r)+V


def V_Yukawa(r, m):
    return -m*(np.exp(-r))/r


def plot_V_eff():
    phi0 = 0  # rad
    v0 = 30.3  # km/s 7181.9  # m/s 10*np.sqrt(G*m/r0)
    m1 = 1.989e30
    m2 = 5.972e24
    m = m1*m2/(m1+m2)
    M = m1+m2
    r0 = 147000000  # (1350+6740)# *1000 # m
    r = np.linspace(0.001, 2*r0, 100)

    l = r0*m*v0

    y = V_r(r, m) + (l * l) / (2 * m * r * r)
    plt.plot(r, y)
    plt.axis([0, 2*r0, -1.5, 1e30])
    plt.axhline(y=0, color='black', linestyle='--')
    #plt.legend(r0)
    plt.show()


if __name__=='__main__':
    plot_V_eff()
