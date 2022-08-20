import numpy as np
from potential import V_r
from aux_tools import polar2cartesian
from plotting import plot

V = V_r

G = 6.6738480 * 10 ** (-11)  # Nm^2/Kg^2
# G = 6.673848 * 10 ** (-2)  # km^3/(kg*s^2)ç
r0 = 3000000 #(1350+6740)# *1000 # m
phi0 = 0  # rad
v0 = 1000  # km/s 7181.9  # m/s 10*np.sqrt(G*m/r0)
m1 = 1e26  # 2300  # kg Satellite
m2 = 1  # 1e26  # 5.972*10**24  # kg Earth
m = m1*m2/(m1+m2)
E = V(r0, G*m) + 1/2*m*v0*v0

l = r0*m*v0


def r(phi):
    p = l*l/(G*m)
    e = np.sqrt(1+2*E*l*l/(G*m*G*m))
    return p/(1+e*np.cos(phi-phi0))


phi = np.linspace(0, 2*np.pi, 1000)

rx, ry, rz = polar2cartesian(r(phi), phi, 0)
a = r(0)
plot(rx, ry, rz)