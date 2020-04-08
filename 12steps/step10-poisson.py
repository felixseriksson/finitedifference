# d2p/dx2 + d2p/dy2 = b
# <=> del^2p = b
# discretize using 2nd order CD
# (pi+1j^n - 2pij^n + pi-1j^n)/dx^2 + (pij+1^n - 2pij^n + pij-1^n)/dy2 = bij^n
# we isolate centrepoint pij^n:
# pij^n = (dy2(pi+1j^n + pi-1j^n) + dx2(pij+1^n + pij-1^n) - bij^n*dx2*dxy)/(2(dx2 + dy2))
# Domain (0, 2) in x and (0, 1) in y
# IC:
# p = 0 everywhere
# BC:
# p = 0 at x = 0 and x = 2, p = 0 at y = 0 and y = 1
# source term: (two spikes)
# bij = 100 at i = nx/4 and j = ny/4
# bij = -100 at i = 3nx/4 and j = 3ny/4
# bij = 0 everywhere else


import numpy
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot2D(x, y, p):
    """ takes as inputs a 1d array of x-values, a 1d array of y values and a 2d array of corresponding pressure values to plot in 3d space"""
    fig = plt.figure(figsize=(11, 7), dpi=100)
    ax = fig.gca(projection='3d')
    X, Y = numpy.meshgrid(x, y)
    ax.plot_surface(X, Y, p[:], rstride=1, cstride=1,
                    cmap="viridis", linewidth=0, antialiased=False)
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 1)
    ax.view_init(30, 225)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')

nx = 50
ny = 50
nt = 100
xmin = 0
xmax = 2
ymin = 0
ymax = 1
dx = (xmax - xmin) / (nx - 1)
dy = (ymax - ymin) / (ny - 1)

p  = numpy.zeros((ny, nx))
pd = numpy.zeros((ny, nx))
b  = numpy.zeros((ny, nx))
x  = numpy.linspace(xmin, xmax, nx)
y  = numpy.linspace(xmin, xmax, ny)

# Source
b[int(ny / 4), int(nx / 4)]  = 100
b[int(3 * ny / 4), int(3 * nx / 4)] = -100

for it in range(nt):

    pd = p.copy()

    p[1:-1,1:-1] = (((pd[1:-1, 2:] + pd[1:-1, :-2]) * dy**2 +
                    (pd[2:, 1:-1] + pd[:-2, 1:-1]) * dx**2 -
                    b[1:-1, 1:-1] * dx**2 * dy**2) / 
                    (2 * (dx**2 + dy**2)))

    p[0, :] = 0
    p[ny-1, :] = 0
    p[:, 0] = 0
    p[:, nx-1] = 0
    if it % 20 == 0:
        plot2D(x, y, p)

plot2D(x, y, p)