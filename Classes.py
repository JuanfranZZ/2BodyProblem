# Inspired in Vallado, Fundamentals in Astrodynamics and applications 4th edition
# Assumptions:
# 1-The mass of the satellite is negligible compared to that
#   of the attracting body. Reasonable for artificial satellites
#   in the foreseeable future.
# 2-The coordinate system chosen is inertial.
# 3-Spherical bodies
# 4-No other forces act on the system except for gravitational
#   forces that act along a line joinning the centers of the two
#   bodies.

import numpy as np
from scipy.integrate import solve_ivp


class Body:
    def __init__(self, name, mass, **kwargs):
        self.name = name
        self.mass = mass
        try:
            if "pos" in kwargs:
                self.set_position(kwargs['pos'])
            else:
                self.pos = None
            if "vel" in kwargs:
                self.set_velocity(kwargs['vel'])
            else:
                self.vel = None
        except Exception as e:
            print(e)
            pass

    def set_position(self, position):
        self.pos = np.array(position)

    def set_velocity(self, velocity):
        self.vel = np.array(velocity)


def body_decoder(obj):
    if "ecc" in obj:
        e = obj["ecc"]
        vel = obj["vel"]*np.sqrt(2/(1-e)-1) #perihelio vel
    else:
        vel = obj["vel"]
    return Body(obj["name"], obj["mass"], pos=obj["pos"], vel=[0, vel, 0])


def unidades_decoder(obj):
    c = {}
    for key, value in obj.items():
        c[str(key)] = value
    return c


def constant_decoder(obj):
    c = {}
    for key, value in obj.items():
        c[str(key)] = value
    return c


class BodySystem:
    def __init__(self, G=None):
        if G:
            self.G = G
        else:
            # self.G = 6.673848*10**(-11) # Nm2/kg2
            self.G = 6.673848*10**(-20)  # km^3/(kg*s^2)
        self.bodies = []

    def add_body(self, body):
        self.bodies.append(body)

    def remove_body(self, body):
        self.bodies.remove(body)


class TwoBodySystem(BodySystem):
    def __init__(self, G=None, body1=None, body2=None):
        super().__init__(G)
        self.bodies = [body1, body2]
        self.r = self.bodies[1].pos - self.bodies[0].pos
        self.r_distance = np.linalg.norm(self.r)
        self.v = self.bodies[1].vel - self.bodies[0].vel
        self.vr_module = np.linalg.norm(np.dot(self.r, self.v) / self.r_distance)

        self.reduced_mass = np.prod([x.mass for x in self.bodies])/np.sum(x.mass for x in self.bodies)  # reduced mass from wikipedia
        self.mu = self.G * np.sum(x.mass for x in self.bodies)
        self.spec_angular_momemtum = self._specific_angular_momentum()
        self.spec_angular_momemtum_module = np.linalg.norm(self.spec_angular_momemtum)
        self.spec_mec_energy = self._specific_mechanical_energy()
        self.eccentricity = self._orbit_eccentricity()
        self._tolerance = 0.00001

    def set_tolerance(self, tol):
        self._tolerance = tol

    def get_tolerance(self):
        return self._tolerance

    def calculate_orbit(self, theta0, thetafin):
        # take initial conditions
        r0 = self.r_distance
        drdt0 = self.vr_module  # radial velocity projection
        h = self.spec_angular_momemtum_module
        u0 = [1/r0, drdt0/(-h)]

        # theta_span = np.linspace(theta0, thetafin, int(1/self.get_tolerance()))
        theta_span = (theta0, thetafin)

        # define dynamic equations
        def dynamic(theta, u):
            y = np.zeros((2, 1)).reshape((2,))
            y[0] = u[1]
            dduddtheta = -u[0] + self.mu/(h*h)
            y[1] = dduddtheta
            return y

        sol = solve_ivp(dynamic, theta_span, u0, method='RK45', max_step=self.get_tolerance())
        sol_r = 1/sol.y[0]
        sol_theta = sol.t
        return sol_r, sol_theta

    def _orbit_eccentricity(self):
        return np.sqrt(1 + 2*self.spec_angular_momemtum_module**2*self.spec_mec_energy/(self.mu**2))
        # return np.linalg.norm(np.cross(self.v, self.spec_angular_momemtum)/self.mu - self.r/np.linalg.norm(self.r))

    def _specific_mechanical_energy(self):
        return self.vr_module ** 2 / 2 + self.V_ef()

    def V(self):
        return - self.mu/self.r_distance

    def V_ef(self):
        return 1/2*(self.spec_angular_momemtum_module/self.r_distance)**2 + self.V()

    def _specific_angular_momentum(self):
        return np.cross(self.r, self.v)
