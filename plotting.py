from matplotlib import pyplot as plt
import numpy as np
from aux_tools import polar2cartesian
from dynamics import theoretical_orbit


def plot(rx, ry, rz, title="", excen=None, mu=None, h=None, theta_fin=None):
    # 3d plot
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)#, projection='3d')

    # ax.view_init(elev=90, azim=-90)
    ax.set_aspect('equal')

    # plot trajectory and starting point

    ax.plot(rx[0], ry[0], 'ro')#, rz[0], 'ro')
    ax.plot(0, 0, 'bo', markersize=14)# 0, 'bo', markersize=14)

    ax.plot(rx, ry, 'k-')#, rz, 'k-')

    # Theoretical orbit
    if excen and mu and h and theta_fin:
        r_theo = theoretical_orbit(h, mu, excen, np.linspace(0, theta_fin, rx.size))
        theta_theo = np.linspace(0, theta_fin, rx.size)
        rx_theo, ry_theo, rz_theo = polar2cartesian(r_theo, theta_theo, 0)
        if rx_theo[0] - rx[0] > rx_theo[int(rx_theo.size/2)] - rx[int(rx.size/2)]:
            theta_theo = np.linspace(np.pi, theta_fin + np.pi, rx.size)
            rx_theo, ry_theo, rz_theo = polar2cartesian(r_theo, theta_theo, 0)
        ax.plot(rx_theo, ry_theo, rz_theo, 'g--')

    plt.legend(['Posicion inicial', 'Punto focal', 'Numérica', 'Teórica'])
    plt.title(title)

    plt.show()

    fig2 = plt.figure(figsize=(10, 10))
    error = np.sqrt((rx-rx_theo)**2 + (ry-ry_theo)**2)
    plt.plot(error)
    plt.show()
