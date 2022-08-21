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
    return Body(obj["name"], obj["mass"], pos=obj["pos"], vel=obj["vel"])


class BodySystem:
    def __init__(self):
        # self.G = 6.673848*10**(-11) # Nm2/kg2
        self.G = 6.673848*10**(-20)  # km^3/(kg*s^2)
        self.bodies = []

    def add_body(self, body):
        self.bodies.append(body)

    def remove_body(self, body):
        self.bodies.remove(body)


class TwoBodySystem(BodySystem):
    def __init__(self, body1, body2):
        super().__init__()
        self.bodies = [body1, body2]
        self.r = self.bodies[0].pos - self.bodies[1].pos
        self.r_distance = np.linalg.norm(self.r)
        self.v = self.bodies[0].vel - self.bodies[1].vel
        self.v_module = np.linalg.norm(self.v)

        self.mu = self.G * np.sum(x.mass for x in self.bodies)
        self.eccentricity = self._orbit_eccentricity()
        self.spec_mec_energy = self._specific_mechanical_energy()
        self.spec_angular_momemtum = self._specific_angular_momentum()
        self._tolerance = 0.00001

    def set_tolerance(self, tol):
        self._tolerance = tol

    def get_tolerance(self):
        return self._tolerance

    def calculate_orbit(self, theta0, thetafin):
        # take initial conditions
        r0 = self.r_distance
        drdt0 = np.dot(self.v, self.r) / self.r_distance  # radial velocity projection
        h = self.spec_angular_momemtum
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
        return (np.dot((self.v_module**2 - self.mu/self.r_distance), self.r) - np.dot(np.dot(self.r, self.v), self.v)) / self.mu

    def _specific_mechanical_energy(self):
        return self.v_module**2/2 - self.mu/self.r_distance

    def _specific_angular_momentum(self):
        return np.linalg.norm(np.cross(self.r, self.v))
