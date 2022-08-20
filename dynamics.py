# MOVEMENT EQUATIONS:
import numpy as np


def odes1(t, x, l, m1, m2, E, V, G, tfin):
    # m1 = 2300  # kg
    r = x[0]
    phi = x[1]
    M = m1+m2
    m = m1*m2/M
    mu = G*M

    # define system of equation
    dphiDt = l / (m * r * r)
    aux = 2 / m * (E - V(r, G*m*M) - (l * l) / (2 * m * r * r))
    drDt = -np.sqrt(aux)
    # print('sqrt =', 2 / m * (E - V(r, m) - (l * l) / (2 * m * r * r)))
    e = np.sqrt(1 + E * l * l / (mu * mu * m * m * m))

    p = True

    if p:
        print('t=', tfin - t)
        print('r=', r)
        print('phi=', phi)
        print('e=', e)
        print('l=', l)
        print('E=', E)
        print('V=', V(r, mu*m))
        print('Veff=', V(r, mu*m) + (l * l) / (2 * m * r * r))
        print('aux=', aux)
        print('=============')

    return [drDt, dphiDt]


def odes2(t, x, l, m1, m2, E, V, G, tfin):
    # m1 = 2300  # kg
    m = m1*m2/(m1+m2)
    r = x[0]
    phi = x[1]
    print('t=', t)
    #print('x=', x)
    #print('l=', l)
    #print('m=', m)
    print('E=', E)
    print('V=', V(r, G*m))
    print('V_eff=', V(r, G*m) + (l * l) / (2 * m * r * r))
    print('=============')

    # define system of equation
    dphiDr = l / (r*r) * 1 / (np.sqrt(2 * m * (E - V(r, m) - (l * l) )))
    drDt = np.sqrt(2 / m * (E - V(r, G*m) - (l * l) / (2 * m * r * r)))
    #print('sqrt =', 2 / m * (E - V(r, m) - (l * l) / (2 * m * r * r)))

    return [drDt, dphiDr]


def odes3(theta, x, h, m1, m2, E, V, G):
    # m1 = 2300  # kg
    u = x[0]
    duDtheta = x[1]
    M = m1+m2
    m = m1*m2/M
    mu = G*M
    dthetaDt = h * (u*u)
    h = dthetaDt/(u*u)

    # define system of equation
    d2uDtheta2 = mu/(h*h) - u
    e = np.sqrt(1 + E * h * h / (mu * mu * m))

    p = True

    if p:
        print('r=', 1/u)
        print('theta=', theta)
        print('e=', e)
        print('l=', h*m)
        print('E=', E)
        print('V=', V(1/u, mu*m))
        print('Veff=', V(1/u, mu*m) + (h*h*m*m) / (2 * m / (u * u)))
        print('=============')

    return [duDtheta, d2uDtheta2]


def theoretical_orbit(h, mu, e, theta):
    r = (h*h/mu)/(1+e*(e*np.cos(theta)))
    return r