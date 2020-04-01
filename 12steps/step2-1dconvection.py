# a k a inviscid burgers'
# can generate (something similar to) shocks
# discontinuity from smooth starting conditions

# discretize: FD in time, BD in space
# ui^n+1 = ui^n - ui^n(dt/dx)(ui^n-ui-1^n) (^n means time index n)
# Initial conds. u = 2 between x = 0.5 and x = 1
# u = 1 elsewhere on [0,2]

import numpy as np
from matplotlib import pyplot as plt
import time
import sys

def steps(nx, sigma):
    dx = 2/(nx-1)
    nt = 20
    c = 1
    dt = sigma * dx

    return nx, dx, nt, dt, c


nx, dx, nt, dt, c = steps(41, 0.5)

u = np.ones(nx)
u[int(0.5/dx):int(1/dx+1)] = 2
print(u)

plt.plot(np.linspace(0, 2, nx), u)
plt.show()

un = np.ones(nx)

for n in range(nt):
    un = u.copy()
    for i in range(1, nx):
        u[i] = un[i] - un[i]* (dt/dx)* (un[i]- un[i-1])
    if n % 10 == 0:
        plt.plot(np.linspace(0, 2, nx), u)
        plt.show()

plt.plot(np.linspace(0, 2, nx), u)
plt.show()