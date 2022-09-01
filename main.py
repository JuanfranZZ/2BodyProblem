import json
import numpy as np
from matplotlib import pyplot as plt
from plotting import plot, polar2cartesian
from Classes import TwoBodySystem, body_decoder, constant_decoder

# leer fichero de planetas
planets_file = r"bodies.json"
with open(planets_file) as file:
    cuerpos = json.loads(file.read()) # https://starlust.org/the-planets-in-order-from-the-sun/

# crear objetos de planetas
# for s in planets:
    # exec(f"{s} = body_decoder(planets['{s}'])")

Tierra = body_decoder(cuerpos['planetas']['Tierra'])
Tierra.set_velocity(Tierra.vel*0.001)
Tierra.set_position(2*Tierra.pos)
Tierra2 = body_decoder(cuerpos['planetas']['Tierra'])
Tierra2.set_position([0,0,0])
Tierra2.set_velocity([0,0,0])
Jupiter = body_decoder(cuerpos['planetas']['Jupiter'])
Venus = body_decoder(cuerpos['planetas']['Venus'])
Sol = body_decoder(cuerpos['Sol'])
constants = constant_decoder(cuerpos["constantes"])
G = constants["G"]["value"]

# Jupiter Earth system
theta0 = 0
thetafin = 2*np.pi
twoBodySystem = TwoBodySystem(G=G, body1=Tierra2, body2=Tierra)
r0 = twoBodySystem.r_distance
vr_0 = twoBodySystem.vr_module
h0 = twoBodySystem.spec_angular_momemtum_module
mu = twoBodySystem.mu
E = twoBodySystem.spec_mec_energy
excen = twoBodySystem.eccentricity

print("r0:", r0)
print("vr_0:", vr_0)
print("h0:", h0)
print("mu:", mu)
print("E:", E)
print("excentricidad:", excen)

# Calcular
orbit = twoBodySystem.calculate_orbit(theta0, thetafin)
r = orbit[0]
theta = orbit[1]

# pintar resultados
# theta = np.linspace(theta0, thetafin, np.size(r))
plt.plot(r)
plt.show()
[rx, ry, rz] = polar2cartesian(r, theta, 0)
title = "r0="+str(r0)+"; v0="+str(vr_0)
plot(rx, ry, rz, title=title, excen=excen, h=h0, mu=mu, theta_fin=thetafin)
