from plotting import plot
from aux_tools import polar2cartesian
from potential import V_r, V_no_m,V_eff, V_Yukawa

import numpy as np
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

V = V_Yukawa


# MOVEMENT EQUATIONS:
def odes(t, x, l, m, E):
    # m1 = 2300  # kg
    r = x[0]
    phi = x[1]
    print('t=', t)
    print('x=', x)
    print('l=', l)
    print('m=', m)
    print('E=', E)
    print('V=', V(r, m))
    print('=============')

    # define system of equation
    dphiDt = l / (m * r * r)
    drDt = np.sqrt(2 / m * (E - V(r, m) - (l * l) / (2 * m * r * r)))
    print('sqrt =', 2 / m * (E - V(r, m) - (l * l) / (2 * m * r * r)))

    return [drDt, dphiDt]


G = 6.6738480 * 10 ** (-11)  # Nm^2/Kg^2

# Initial conditions

r0 = (1350+6740)*1000 # m
phi0 = 0  # rad
v0 = 7181.9  # m/s 10*np.sqrt(G*m/r0)

# Masses

m1 = 2300  # kg Satellite
m2 = 5.972*10**24  # kg Earth
m = m1*m2/(m1+m2)
print('m=', m)

# Angular moment

l = r0*m*v0

# Tukawa example

m = 1
l = np.sqrt(0.5)
r0 = 0.4
v0 = l/(m*r0)

print('V=', V(r0, m))
E0 = V(r0, m) + 1/2*m*v0*v0
print('E0=', E0)
print('K0=', 1/2*m*v0*v0)
print('---------------------------')

tini = 0
tfin = 10
tolerance = 0.0001
op = 2

t = np.arange(tini, tfin, tolerance)  # for ideint
t_span = (tini, tfin)  # for solve_ivp

if op == 1:  # odeint
    sol = odeint(odes, (r0, phi0), t, args=(l, m, E0), tfirst=True, full_output=1)
    r = sol[0][:, 0]
    phi = sol[0][:, 1]
elif op == 2:  # RK45
    sol = solve_ivp(odes, t_span, (r0, phi0), args=(l, m, E0), method='RK45', max_step=tolerance)
    r = sol.y[0]
    phi = sol.y[1]
elif op == 3:  # LSODA
    sol = solve_ivp(odes, t_span, (r0, phi0), args=(l, m, E0), method='LSODA', max_step=tolerance)
    r = sol.y[0]
    phi = sol.y[1]



[rx, ry, rz] = polar2cartesian(r, phi, 0)
plot(rx, ry, rz)