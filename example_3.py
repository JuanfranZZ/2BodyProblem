import numpy as np

G = 6.67430e-20  # km**3/(kg * s**2)
m_1 = m_2 = 1.0e26  # kg

R_1_0 = np.array((0, 0, 0))  # km
R_2_0 = np.array((3000, 0, 0))  # km
dotR_1_0 = np.array((10, 20, 30))  # km/s
dotR_2_0 = np.array((0, 40, 0))  # km/s

y_0 = np.hstack((R_1_0, R_2_0, dotR_1_0, dotR_2_0))

X_1 = y_0[0]
Y_1 = y_0[1]
Z_1 = y_0[2]
X_2 = y_0[3]
Y_2 = y_0[4]
Z_2 = y_0[5]

r = np.sqrt((X_2 - X_1) ** 2 + (Y_2 - Y_1) ** 2 + (Z_2 - Z_1) ** 2)

ddotX_1 = G * m_2 * (X_2 - X_1) / r ** 3
ddotY_1 = G * m_2 * (Y_2 - Y_1) / r ** 3
ddotZ_1 = G * m_2 * (Z_2 - Z_1) / r ** 3
ddotX_2 = -G * m_1 * (X_2 - X_1) / r ** 3
ddotY_2 = -G * m_1 * (Y_2 - Y_1) / r ** 3
ddotZ_2 = -G * m_1 * (Z_2 - Z_1) / r ** 3

X_1 = y_0[0]
Y_1 = y_0[1]
Z_1 = y_0[2]
X_2 = y_0[3]
Y_2 = y_0[4]
Z_2 = y_0[5]

r = np.sqrt((X_2 - X_1) ** 2 + (Y_2 - Y_1) ** 2 + (Z_2 - Z_1) ** 2)

ddotX_1 = G * m_2 * (X_2 - X_1) / r ** 3
ddotY_1 = G * m_2 * (Y_2 - Y_1) / r ** 3
ddotZ_1 = G * m_2 * (Z_2 - Z_1) / r ** 3
ddotX_2 = -G * m_1 * (X_2 - X_1) / r ** 3
ddotY_2 = -G * m_1 * (Y_2 - Y_1) / r ** 3
ddotZ_2 = -G * m_1 * (Z_2 - Z_1) / r ** 3

R_1 = y_0[:3]
R_2 = y_0[3:6]

r = np.sqrt(np.sum(np.square(R_2 - R_1)))
ddot = G * (R_2 - R_1) / r ** 3
ddotR_1_0 = m_2 * ddot
ddotR_2_0 = -m_1 * ddot

Delta_t = 1  # s
dotR_1_1 = ddotR_1_0 * Delta_t + dotR_1_0
dotR_2_1 = ddotR_2_0 * Delta_t + dotR_2_0

R_1_1 = dotR_1_0 * Delta_t + R_1_0
R_2_1 = dotR_2_0 * Delta_t + R_2_0