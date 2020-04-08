# FD in time, BD in space (as with 2d burgers), CD for pressure (usual difference for laplacian)
# discretized:
# LHS = (uij^n+1 - uij^n)/dt + uij^n((uij^n - ui-1j^n)/dx) + vij^n((uij^n - uij-1^n)/dy) = 
# RHS = -1/rho((pi+1j^n - pi-1j^n)/(2dx)) + nu((ui+1j^n - 2uij^n + ui-1j^n)/dx2 + (uij+1^n - 2uij^n + uij-1^n)/dy2)
# same for v also...
# poisson:
# LHS = (pi+1j^n - pij^n + pi-1j^n)/dx2 + (pij+1^n - pij^n + pij-1^n)/dy2 = 
# RHS = -rho(((ui+1j^n - ui-1j^n)/(2dx))^2 + 2((uij+1^n - uij-1^n)/(2dy))((vi+1j^n - vi-1j^n)/(2dx)) + ((vij+1^n - vij-1^n)/(2dy))^2) + rho/dt((ui+1j^n - ui-1j^n)/(2dx) + (vij+1^n - vij-1^n)/(2dy))
# transpose by isolating, uij^n+1, vij^n+1 and pij^n
# ICs:
# u, v, p = 0 everywehre
# BCs:
# u = 1 at y = 2
# u, v = 0 at x = 0, 2 and y = 0 (basically no-slip everywhere except at y = 2)
# p = 0 at y = 2
# dp/dy = 0 at y = 0 and dp/dx = 0 at x = 0, 2
# we will have actual time steps - nt - and artificial time steps - nit - for iteration of pressure inside outer loop

import numpy
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

nx = 41
ny = 41
nt = 500
nit = 50
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)
X, Y = numpy.meshgrid(x, y)

rho = 1
nu = 0.1
dt = 0.001

u = numpy.zeros((ny, nx))
v = numpy.zeros((ny, nx))
p = numpy.zeros((ny, nx))
b = numpy.zeros((ny, nx))

def build_up_b(b, rho, dt, u, v, dx, dy):
    b[1:-1, 1:-1] = (rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) + (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) - 
    ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 - 2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) * 
    (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)) - ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2))

    return b

def pressure_poisson(p, dx, dy, b):
    pn = numpy.empty_like(p)
    pn = p.copy()
    
    for q in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 + (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) / (2 * (dx**2 + dy**2)) - 
        dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1,1:-1])

        p[:, -1] = p[:, -2] # dp/dx = 0 at x = 2
        p[0, :] = p[1, :]   # dp/dy = 0 at y = 0
        p[:, 0] = p[:, 1]   # dp/dx = 0 at x = 0
        p[-1, :] = 0        # p = 0 at y = 2
        
    return p

def cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu):
    un = numpy.empty_like(u)
    vn = numpy.empty_like(v)
    b = numpy.zeros((ny, nx))
    
    for n in range(nt):
        un = u.copy()
        vn = v.copy()
        
        b = build_up_b(b, rho, dt, u, v, dx, dy)
        p = pressure_poisson(p, dx, dy, b)
        
        u[1:-1, 1:-1] = (un[1:-1, 1:-1]-
                         un[1:-1, 1:-1] * dt / dx *
                        (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                         vn[1:-1, 1:-1] * dt / dy *
                        (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                         dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                         nu * (dt / dx**2 *
                        (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                         dt / dy**2 *
                        (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))

        v[1:-1,1:-1] = (vn[1:-1, 1:-1] -
                        un[1:-1, 1:-1] * dt / dx *
                       (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                        vn[1:-1, 1:-1] * dt / dy *
                       (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                        dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) +
                        nu * (dt / dx**2 *
                       (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                        dt / dy**2 *
                       (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

        u[0, :]  = 0
        u[:, 0]  = 0
        u[:, -1] = 0
        u[-1, :] = 1    # set velocity on cavity lid equal to 1
        v[0, :]  = 0
        v[-1, :] = 0
        v[:, 0]  = 0
        v[:, -1] = 0
        
        
    return u, v, p

# main code 

u = numpy.zeros((ny, nx))
v = numpy.zeros((ny, nx))
p = numpy.zeros((ny, nx))
b = numpy.zeros((ny, nx))
nt = 100
u, v, p = cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu)

fig = plt.figure(figsize=(11,7), dpi=100)
# plotting the pressure field as a contour
plt.contourf(X, Y, p, alpha=0.5, cmap="viridis")  
plt.colorbar()
# plotting the pressure field outlines
plt.contour(X, Y, p, cmap="viridis")  
# plotting velocity field
plt.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2]) 
plt.xlabel('X')
plt.ylabel('Y')

# Notes:
# We can test with different values of nt to see the system stabilise over time
# We can also plot using "plt.streamplot(X, Y, u, v)" instead of quiver of every second vector