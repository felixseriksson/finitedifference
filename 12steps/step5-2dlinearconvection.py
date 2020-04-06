# extend from 1d to 2d: apply multivariable calculus rules, build 2d grid and discretize
# by keeping x, then y constant
# uij = u(xi, yj)
# ui+1j+1 = uij + (dx*d/dx + dy*d/dy)(uij) + hot (simple taylor series expansion)
# can do FD, BD and CD in x, y and time
# du/dt + cdu/dx + cdu/dy = 0
# FD in time, BD in x and y:
# uij^n+1 = uij^n - cdt/dx(uij^n - ui-1j^n) - cdt/dy(uij^n-uij-1^n)
# ICs/BCs:
# square wave but in 3d, so cube (?) wave? hypersquare wave ?


from mpl_toolkits.mplot3d import Axes3D

import numpy
from matplotlib import pyplot as plt

###variable declarations
nx = 81
ny = 81
nt = 100
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .2
dt = sigma * dx

x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)

u = numpy.ones((ny, nx))
un = numpy.ones((ny, nx))

# initial conditions

##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2 

###Plot Initial Condition
##the figsize parameter can be used to produce different sized images

fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = numpy.meshgrid(x, y)
surf = ax.plot_surface(X, Y, u[:], cmap="viridis")
plt.show()

# Test using nested loops
'''
u = numpy.ones((ny, nx))
u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2

for n in range(nt + 1): ##loop across number of time steps
    un = u.copy()
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
'''

# Test using array operations
u = numpy.ones((ny, nx))
u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2

for n in range(nt + 1): ##loop across number of time steps
    un = u.copy()
    u[1:, 1:] = (un[1:, 1:] - (c * dt / dx * (un[1:, 1:] - un[1:, :-1])) - (c * dt / dy * (un[1:, 1:] - un[:-1, 1:])))
    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1
    if n % 20 == 0:
        fig = plt.figure(figsize=(11, 7), dpi=100)
        ax = fig.gca(projection='3d')
        surf2 = ax.plot_surface(X, Y, u[:], cmap="viridis")
        plt.show()

fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
surf2 = ax.plot_surface(X, Y, u[:], cmap="viridis")
plt.show()