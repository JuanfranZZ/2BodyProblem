from matplotlib import pyplot as plt
import numpy as np


def plot(rx, ry, rz):
    # 3d plot
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    ax.view_init(elev=90, azim=-90)

    # plot trajectory and starting point
    ax.plot(rx, ry, rz)
    ax.plot(rx[0], ry[0], rz[0], 'ro')
    ax.plot(0, 0, 0, 'bo', markersize=14)

    plt.legend(['Trajectory', 'Starting Position'])

    plt.show()