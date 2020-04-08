# d2p/dx2 + d2p/dy2 = 0
# <=> del^2p = 0
# discretize using 2nd order CD
# (pi+1j^n - 2pij^n + pi-1j^n)/dx^2 + (pij+1^n - 2pij^n + pij-1^n)/dy2 = 0
# we isolate centrepoint pij^n:
# pij^n = (dy2(pi+1j^n + pi-1j^n) + dx2(pij+1^n + pij-1^n))/(2(dx2 + dy2))
# Domain (0, 2) in x and (0, 1) in y
# IC:
# p = 0 everywhere
# BC:
# p = 0 at x = 0
# p = y at x = 2
# dp/dy = 0 at y = 0 and 1
# Analytical sol: (can be obtained by separating variables)
# p(x, y) = x/4 - 4*infsum over odd n (sinh(npix)cosh(npiy)/(n^2pi^2*sinh(2pin)))
# Notes:
# 1) laplacian is typical of diffusion (isotropic) --> has to be discretized with CD, to be consistent with physics
# 2) 2nd order CD in x and y is most used num. scheme for laplacian, also known as five point difference operator or five point stencil
# 3) iterative method for a steady state (no time), (artificial time variable) --> linear system of eqs. with pentadiagonal coeff. matrix
# equivalent to the "point jacobi method"

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


def laplace2d(p, y, dx, dy, l1norm_target):
    """takes 2d p, 1d y, stepsizes dy and dx as well as a criterion for convergence, iterates until criterion is met"""
    l1norm = 1
    pn = numpy.empty_like(p)

    while l1norm > l1norm_target:
        pn = p.copy()
        p[1:-1, 1:-1] = ((dy**2 * (pn[1:-1, 2:] + pn[1:-1, 0:-2]) + dx ** 2 * (pn[2:, 1:-1] + pn[0:-2, 1:-1])) / (2 * (dx**2 + dy**2)))

        p[:, 0] = 0  # p = 0 @ x = 0
        p[:, -1] = y  # p = y @ x = 2
        p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
        p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = 1
        l1norm = (numpy.sum(numpy.abs(p[:]) - numpy.abs(pn[:])) / numpy.sum(numpy.abs(pn[:])))

    return p


nx = 31
ny = 31
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)

p = numpy.zeros((ny, nx))
x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 1, ny)

# BCs
p[:, 0] = 0  # p = 0 @ x = 0
p[:, -1] = y  # p = y @ x = 2
p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = 1

plot2D(x, y, p)

p = laplace2d(p, y, dx, dy, 1e-4)

plot2D(x, y, p)