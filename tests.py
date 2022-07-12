import numpy as np
from potential import V
from matplotlib import pyplot as plt

m1 = 2300 # kg
m2 = 5.972*10**24 # kg
m = m1*m2/(m1+m2)

r = np.linspace(6370000, 8000000)

plt.plot(r, V(r, m))
plt.show()
