from plotting import plot
from aux_tools import polar2cartesian
from potential import V_adim

import numpy as np
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

V = V_adim


# MOVEMENT EQUATIONS:
def odes(t, x, l, m, E):
    m1 = 2300  # kg
    r = x[0]
    phi = x[1]
    print('t=', t)
    print('x=', x)
    print('l=', l)
    print('m=', m)
    print('E=', E)
    print('=============')

    # define system of equation
    dphiDt = l / (m * r * r)
    drDt = np.sqrt(2 / m1 * (E - V(r, m) - (l * l) / (2 * m * r * r)))
    print('sqrt =', 2 / m1 * (E - V(r, m) - (l * l) / (2 * m * r * r)))

    return [drDt, dphiDt]


G = 6.6738480 * 10 ** (-11)  # Nm^2/Kg^2

# Initial conditions

r0 = (1350+6740)*1000 # m
phi0 = 0  # rad

m1 = 2300  # kg
m2 = 5.972*10**24  # kg
m = m1*m2/(m1+m2)
v0 = 7181.9# m/s 10*np.sqrt(G*m/r0)
print('m=', m)

l = r0*m*v0
print('V=', V(r0, m))
E0 = V(r0, m)+1/2*m*v0*v0
print('E0=', E0)
print('K0=', 1/2*m*v0*v0)
print('---------------------------')

tini = 0
tfin = 3600*24
tolerance = 0.1
op = 1

t = np.arange(tini, tfin, tolerance)  # for ideint
t_span = (tini, tfin)  # for solve_ivp

if op == 1:
    sol = odeint(odes, (r0, phi0), t, args=(l, m, E0), tfirst=True, full_output=1)
    r = sol[0][:, 0]
    phi = sol[0][:, 1]
elif op == 2:
    sol = solve_ivp(odes, t_span, (r0,phi0), args=(l, m, E0), method='RK45', t_eval=t)
    r = sol.y[0]
    phi = sol.y[1]


[rx, ry, rz] = polar2cartesian(r, phi, 0)
plot(rx, ry, rz)