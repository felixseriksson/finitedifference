# FD in time, BD in x and y
# bit harder to discretize in 2d
# because we have 2 diff. eqs. now
# du/dt + u(du/dx) + v(du/dy) = 0
# dv/dt + u(dv/dx) + v(dv/dy) = 0
# discretizing gives:
# uij^n+1 = uij^n - uij^n(dt/dx)(uij^n - ui-1j^n) - vij^n(dt/dy)(uij^n - uij-1^n)
# vij^n+1 = vij^n - uij^n(dt/dx)(vij^n - vi-1j^n) - vij^n(dt/dy)(vij^n - vij-1^n)
# ICs:
# same, 2d hat function
# BCs:
# also same

from mpl_toolkits.mplot3d import Axes3D

import numpy
from matplotlib import pyplot as plt

###variable declarations
nx = 101
ny = 101
nt = 80
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .2
dt = sigma * dx

x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)

u = numpy.ones((ny, nx))
un = numpy.ones((ny, nx))
v = numpy.ones((ny, nx))
vn = numpy.ones((ny, nx))

# initial conditions

##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2
##set hat function I.C. : v(.5<=x<=1 && .5<=y<=1 ) is 2
v[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2

###Plot Initial Condition
##the figsize parameter can be used to produce different sized images

fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = numpy.meshgrid(x, y)
surf = ax.plot_surface(X, Y, u[:], cmap="viridis", rstride=2, cstride=2)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
plt.show()

# Test using nested loops
'''
# NOT UPDATED, WILL PROBABLY NOT WORK AS INTENDED

for n in range(nt + 1): ##loop across number of time steps
    un = u.copy()
    vn.copy()
    row, col = u.shape
    for j in range(1, row):
        for i in range(1, col):
            u[j, i] = (un[j, i] - (c * dt / dx * (un[j, i] - un[j, i - 1])) - (c * dt / dy * (un[j, i] - un[j - 1, i])))
            u[0, :] = 1
            u[-1, :] = 1
            u[:, 0] = 1
            u[:, -1] = 1

fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
surf2 = ax.plot_surface(X, Y, u[:], cmap="viridis")
plt.show()
#'''

# Test using array operations

for n in range(nt + 1): ##loop across number of time steps
    un = u.copy()
    vnn = v.copy()
    u[1:, 1:] = (un[1:, 1:] - (un[1:, 1:] * c * dt / dx * (un[1:, 1:] - un[1:, :-1])) - vn[1:, 1:] * c * dt / dy * (un[1:, 1:] - un[:-1, 1:]))
    v[1:, 1:] = (vn[1:, 1:] - (un[1:, 1:] * c * dt / dx * (vn[1:, 1:] - vn[1:, :-1])) - vn[1:, 1:] * c * dt / dy * (vn[1:, 1:] - vn[:-1, 1:]))
    
    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1
    
    v[0, :] = 1
    v[-1, :] = 1
    v[:, 0] = 1
    v[:, -1] = 1

    if n % 20 == 0:
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
surf2 = ax.plot_surface(X, Y, u, cmap="viridis")
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
plt.show()