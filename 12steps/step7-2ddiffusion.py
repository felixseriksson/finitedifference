# du/dt = nu(d2u/dx2 + d2u/dx2)
# 2nd order CD in space, FD in time
# discretizing
# uij^n+1 = uij^n + nudt/dx^2(ui-1j^n - 2uij^n + ui+1j^n) + nudt/dy^2(uij-1^n - 2uij^n + uij+1^n)
# same ICs, same BCs

from mpl_toolkits.mplot3d import Axes3D

import numpy
from matplotlib import pyplot as plt

###variable declarations
nx = 61
ny = 61
nt = 17
nu = 0.05
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = 0.25
dt = sigma * dx * dy / nu

x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)

u = numpy.ones((ny, nx))
un = numpy.ones((ny, nx))

u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2

# initial conditions

fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = numpy.meshgrid(x, y)
surf = ax.plot_surface(X, Y, u, cmap="viridis")#, rstride=2, cstride=2)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_xlim(0, 2)
ax.set_ylim(0, 2)
ax.set_zlim(1, 2.5)
plt.show()

def diffuse(nt):
    """diffuse over nt timesteps"""
    u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2

    for n in range(nt + 1):
        un = u.copy()
        u[1:-1, 1:-1] = (un[1:-1,1:-1] + nu * dt / dx**2 * (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) + nu * dt / dy**2 * (un[2:,1: -1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1]))

        u[0, :] = 1
        u[-1, :] = 1
        u[:, 0] = 1
        u[:, -1] = 1

        if n % 25 == 0:
            fig = plt.figure(figsize=(11, 7), dpi=100)
            ax = fig.gca(projection='3d')
            X, Y = numpy.meshgrid(x, y)
            _ = ax.plot_surface(X, Y, u, cmap="viridis")
            ax.set_xlabel('$x$')
            ax.set_ylabel('$y$')
            ax.set_xlim(0, 2)
            ax.set_ylim(0, 2)
            ax.set_zlim(1, 2.5) 
            plt.show()

# fig = plt.figure(figsize=(11, 7), dpi=100)
# ax = fig.gca(projection='3d')
# X, Y = numpy.meshgrid(x, y)
# surf2 = ax.plot_surface(X, Y, u, cmap="viridis")
# ax.set_xlabel('$x$')
# ax.set_ylabel('$y$')
# plt.show()

diffuse(141)