# 1D Burgers' equation
# du/dt + u*du/dx = mu*d2u/dx2
# FD in time, BD in first order space, CD in 2nd order space
# ui^n+1 = ui^n - ui^n(dt/dx)(ui^n - ui-1^n) + nu(dt/dx2)(ui+1^n - 2*ui^n + ui-1^n)
# ICs
# u = -2*nu*((dphi/dx)/phi) + 4
# where phi = exp(-x^2/(4nu)) + exp((-(x-2pi)^2)/(4nu))
# BCs
# u(0) = u(2pi) (periodic with T = 2pi)
# analytical sol
# u = -2*nu*((dphi/dx)/phi) + 4
# where phi = exp(-(x-4t)^2/(4nu(t+1))) + exp(-(x-4t-2pi)^2/(4nu(t+1)))

import numpy as np
from matplotlib import pyplot as plt
import sympy

from sympy.utilities.lambdify import lambdify

x, nu, t = sympy.symbols('x nu t')
phi = (sympy.exp(-(x - 4 * t)**2 / (4 * nu * (t + 1))) +
       sympy.exp(-(x - 4 * t - 2 * sympy.pi)**2 / (4 * nu * (t + 1))))
# print(phi)
phiprime = phi.diff(x)
# print(phiprime)

u = -2 * nu * (phiprime/phi) + 4
ufun = lambdify((t, x, nu), u)

# variable decls.
nx = 101
nt = 300
dx = 2 * np.pi / (nx-1)
nu = 0.07
dt = dx * nu

x = np.linspace(0, 2*np.pi, nx)
un = np.empty(nx)
t = 0

u = np.asarray([ufun(t, x0, nu) for x0 in x])
print(u)

plt.figure(figsize=(8, 5), dpi=100)
plt.plot(x, u, marker='.', lw=1, label="Comp.")
plt.xlim([0, 2 * np.pi])
plt.ylim([0, 10])
plt.legend()
plt.show()


for n in range(nt):
    un = u.copy()
    for i in range(1, nx-1):
        u[i] = un[i] - un[i] * dt / dx *(un[i] - un[i-1]) + nu * dt / dx**2 *(un[i+1] - 2 * un[i] + un[i-1])
    u[0] = un[0] - un[0] * dt / dx * (un[0] - un[-2]) + nu * dt / dx**2 *(un[1] - 2 * un[0] + un[-2])
    u[-1] = u[0]
    if n % 25 == 0:
        plt.figure(figsize=(8, 5), dpi=100)
        plt.xlim([0, 2 * np.pi])
        plt.ylim([0, 10])
        plt.plot(x, u, marker='.', lw=1, label="Comp.")
        plt.show()

u_analytical = np.asarray([ufun(nt * dt, xi, nu) for xi in x])

plt.figure(figsize=(8, 5), dpi=100)
plt.plot(x, u, marker='.', lw=1, label="Comp.")
plt.plot(x, u_analytical, label="Analyt.")
plt.xlim([0, 2 * np.pi])
plt.ylim([0, 10])
plt.legend()
plt.show()