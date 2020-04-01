# du/dt = nu*d2u/dx2
# heat equation of u is t
# exact sols known for constant nu
# kan ansätta vågfunktion och få u = uhatt e^ikx*e^-nuk^2t
# exponentiell dämpning i tid av våg i rymd
# (med nu>0)
# diffusion is isotropic => central difference is appropriate!
# Scheme:
# FD in time, CD in space
# ui^n+1 = ui^n + nu(dt/dx^2)(ui+1^n-2ui^n+ui-1^n)
# same ICs and BCs as before

import numpy as np
from matplotlib import pyplot as plt
import time
import sys

'''
def steps(nx, sigma):
    dx = 2/(nx-1)
    nt = 20
    c = 1
    dt = sigma * dx

    return nx, dx, nt, dt, c


nx, dx, nt, dt, c = steps(41, 0.5)
'''

nx = 41
dx = 2 / (nx - 1)
nt = 200    #the number of timesteps we want to calculate
nu = 0.3   #the value of viscosity
sigma = 0.2 #sigma is a parameter, we'll learn more about it later
dt = sigma * dx**2 / nu #dt is defined using sigma ... more later!


u = np.ones(nx)
u[int(0.5/dx):int(1/dx+1)] = 2
print(u)

plt.plot(np.linspace(0, 2, nx), u)
plt.show()

un = np.ones(nx)

for n in range(nt):
    un = u.copy()
    for i in range(1, nx-1):
        u[i] = un[i] + nu * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1])
    if n % 10 == 0:
        plt.plot(np.linspace(0, 2, nx), u)
        plt.show()

plt.plot(np.linspace(0, 2, nx), u)
plt.show()