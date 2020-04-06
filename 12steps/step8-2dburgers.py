# du/dt + u*du/dx + v*du/dy = nu(d2u/dx2 + d2u/dy2)
# dv/dt + u*dv/dx + v*dv/dy = nu(d2v/dx2 + d2v/dy2)
#
# again, FD in time, BD in x and y and CD in 2nd order x and y
# uij^n+1 = uij^n - uij^n(dt/dx)(uij^n - ui-1j^n) - vij^n(dt/dy)(uij^n - uij-1^n) + nu*dt((ui+1j^n - 2uij^n + ui-1j^n)/dx^2 + (uij+1^n - 2uij^n + uij-1^n)/dy^2)
# for v: same model, but v replaces u in some places
# oh also: same ICs and BCs as before

from mpl_toolkits.mplot3d import Axes3D

import numpy
from matplotlib import pyplot as plt

# variable declarations
nx = 41
ny = 41
nt = 250
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = 0.0009
nu = 0.01
dt = sigma * dx * dy / nu

x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)

u = numpy.ones((ny, nx))
un = numpy.ones((ny, nx))
v = numpy.ones((ny, nx))
vn = numpy.ones((ny, nx))
comb = numpy.ones((ny, nx))

# initial conditions

# set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2
# set hat function I.C. : v(.5<=x<=1 && .5<=y<=1 ) is 2
v[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2

fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = numpy.meshgrid(x, y)
surf = ax.plot_surface(X, Y, u, cmap="viridis", rstride=1, cstride=1)
surf = ax.plot_surface(X, Y, v, cmap="viridis", rstride=1, cstride=1)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
plt.show()

for n in range(nt + 1):  # loop across number of time steps
    un = u.copy()
    vn = v.copy()
    u[1:-1, 1:-1] = (un[1:-1, 1:-1] - dt / dx * un[1:-1, 1:-1] * (un[1:-1, 1:-1] - un[1:-1, 0:-2]) - dt / dy * vn[1:-1, 1:-1] * (un[1:-1, 1:-1] - un[0:-2, 1:-1]) +
                     nu * dt / dx**2 * (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) + nu * dt / dy**2 * (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1]))

    v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - dt / dx * un[1:-1, 1:-1] * (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) - dt / dy * vn[1:-1, 1:-1] * (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) +
                     nu * dt / dx**2 * (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) + nu * dt / dy**2 * (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1]))

    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1

    v[0, :] = 1
    v[-1, :] = 1
    v[:, 0] = 1
    v[:, -1] = 1

    if n % 25 == 0:
        fig = plt.figure(figsize=(11, 7), dpi=100)
        ax = fig.gca(projection='3d')
        X, Y = numpy.meshgrid(x, y)
        surf2 = ax.plot_surface(X, Y, u, cmap="viridis")
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
        plt.show()

fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = numpy.meshgrid(x, y)
ax.plot_surface(X, Y, u, cmap="viridis", rstride=1, cstride=1)
ax.plot_surface(X, Y, v, cmap="viridis", rstride=1, cstride=1)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
plt.show()
