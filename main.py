from plotting import plot
from aux_tools import polar2cartesian
from potential import V_r, V_no_m, V_eff, V_Yukawa
from dynamics import odes1, odes2, odes3

import numpy as np
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from matplotlib import pyplot as plt

V = V_Yukawa

global G

# G = 6.6738480 * 10 ** (-11)  # Nm^2/Kg^2
G = 6.673848 * 10 ** (-20)  # km^3/(kg*s^2)

# Initial conditions

# radio
# ra_sun_earth: 152000 km
# rp_sun_earth: 147000 km

# vel
# earth_Sun perihelio: 30.3km/s
# earth_Sun aphelio: 29.3km/s

r0 = 147000000  # (1350+6740)# *1000 # m
theta0 = 0  # rad
v0 = 30.3  # km/s 7181.9  # m/s 10*np.sqrt(G*m/r0)

# Masses
# Sun: 1.989e30 kg
# Earth: 5.972e24 kg

m1 = 1.989e30
m2 = 5.972e24


# Tukawa example

if V == V_Yukawa:
    m1 = 10000
    m2 = 1
    l = np.sqrt(0.5)
    r0 = 100000000
    m = m1 * m2 / (m1 + m2)
    v0 = l / (m * r0)

m = m1 * m2 / (m1 + m2)
M = m1 + m2
mu = G * M
print('m=', m)

# Angular moment

l0 = r0 * m * v0
K0 = 1 / 2 * m * v0 * v0
E0 = V(r0, mu * m) + K0

# Excentricity

e = np.sqrt(1 + E0 * l0 * l0 / (mu * mu * m * m * m))

print('V=', V(r0, mu * m))
print('E0=', E0)
print('e=', e)
print('K0=', 1 / 2 * m * v0 * v0)
print('l0=', l0)
print('Veff=', V(r0, m * mu) + (l0 * l0) / (2 * m * r0 * r0))
print('---------------------------')

tini = 0
tfin = 60 * 60 * 24 * 365 / 4  # 31536000 s -> 1year

thetafin = np.pi/8

tolerance = 0.0001
op = 2
odes = odes3

#t = np.arange(tini, tfin, tolerance)  # for ideint
#t_span = (tini, tfin)  # for solve_ivp

theta = np.arange(theta0, thetafin, tolerance)
theta_span = (theta0, thetafin)

u0 = 1 / r0
h0 = r0 * v0
duDtheta0 = - np.linalg.norm(v0) / h0

if op == 1:  # odeint
    sol = odeint(odes, (r0, theta0), t, args=(l0, m1, m2, E0, V, G, tfin), tfirst=True, full_output=1)
    r = sol[0][:, 0]
    phi = sol[0][:, 1]
elif op == 2:  # RK45
    sol = solve_ivp(odes, theta_span, (u0, duDtheta0), args=(h0, m1, m2, E0, V, G), method='RK45',
                    max_step=tolerance)
    r = 1/sol.y[0]
    v = 1/sol.y[1]
elif op == 3:  # LSODA
    sol = solve_ivp(odes, t_span, (r0, phi0), args=(l0, m1, m2, E0, V, G, tfin), method='LSODA', max_step=tolerance)
    r = sol.y[0]
    phi = sol.y[1]

theta = np.linspace(theta0, thetafin, np.size(r))
plt.plot(r)
plt.show()
[rx, ry, rz] = polar2cartesian(r, theta, 0)
title = "r0="+str(r0)+"; v0="+str(v0)
plot(rx, ry, rz, title=title, excen=0.7, h=h0, mu=mu, theta_fin=thetafin)
