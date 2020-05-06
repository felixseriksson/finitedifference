import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from grid import generategrid

def pressure_poisson_periodic(p, dx, dy):
    pn = np.empty_like(p)
    
    for _ in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 +
                          (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) /
                         (2 * (dx**2 + dy**2)) -
                         dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, 1:-1])

        # Periodic BC Pressure @ x = 2
        p[1:-1, -1] = (((pn[1:-1, 0] + pn[1:-1, -2])* dy**2 +
                        (pn[2:, -1] + pn[0:-2, -1]) * dx**2) /
                       (2 * (dx**2 + dy**2)) -
                       dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, -1])

        # Periodic BC Pressure @ x = 0
        p[1:-1, 0] = (((pn[1:-1, 1] + pn[1:-1, -1])* dy**2 +
                       (pn[2:, 0] + pn[0:-2, 0]) * dx**2) /
                      (2 * (dx**2 + dy**2)) -
                      dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, 0])
        
        # Wall boundary conditions, pressure
        p[-1, :] =p[-2, :]  # dp/dy = 0 at y = 2
        p[0, :] = p[1, :]  # dp/dy = 0 at y = 0
    
    return p


# variable declarations
nx = 41
ny = 41
nt = 500
nit = 50 
c = 1
displayfreq = 50

# make grid
gridsizex = 16
gridsizey = 4

dx, dy, X, Y = generategrid(gridsizex, gridsizey, nx, ny)

##physical variables
rho = 1
nu = .1
F = 1
dt = .01

#initial conditions
u = np.zeros((ny, nx))
un = np.zeros((ny, nx))

v = np.zeros((ny, nx))
vn = np.zeros((ny, nx))

p = np.ones((ny, nx))
pn = np.ones((ny, nx))

b = np.zeros((ny, nx))

udiff = 1
stepcount = 0

for timestep in range(nt):
    un = u.copy()
    vn = v.copy()

    #b = build_up_b(rho, dt, dx, dy, u, v)
    p = pressure_poisson_periodic(p, dx, dy)

    u[1:-1, 1:-1] = (un[1:-1, 1:-1] -
                    un[1:-1, 1:-1] * dt / dx * 
                    (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                    vn[1:-1, 1:-1] * dt / dy * 
                    (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                    dt / (2 * rho * dx) * 
                    (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                    nu * (dt / dx**2 * 
                    (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                    dt / dy**2 * 
                    (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])) + 
                    F * dt)

    v[1:-1, 1:-1] = (vn[1:-1, 1:-1] -
                    un[1:-1, 1:-1] * dt / dx * 
                    (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                    vn[1:-1, 1:-1] * dt / dy * 
                    (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                    dt / (2 * rho * dy) * 
                    (p[2:, 1:-1] - p[0:-2, 1:-1]) +
                    nu * (dt / dx**2 *
                    (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                    dt / dy**2 * 
                    (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

    # Periodic BC u @ x = 2     
    u[1:-1, -1] = (un[1:-1, -1] - un[1:-1, -1] * dt / dx * 
                (un[1:-1, -1] - un[1:-1, -2]) -
                vn[1:-1, -1] * dt / dy * 
                (un[1:-1, -1] - un[0:-2, -1]) -
                dt / (2 * rho * dx) *
                (p[1:-1, 0] - p[1:-1, -2]) + 
                nu * (dt / dx**2 * 
                (un[1:-1, 0] - 2 * un[1:-1,-1] + un[1:-1, -2]) +
                dt / dy**2 * 
                (un[2:, -1] - 2 * un[1:-1, -1] + un[0:-2, -1])) + F * dt)

    # Periodic BC u @ x = 0
    u[1:-1, 0] = (un[1:-1, 0] - un[1:-1, 0] * dt / dx *
                (un[1:-1, 0] - un[1:-1, -1]) -
                vn[1:-1, 0] * dt / dy * 
                (un[1:-1, 0] - un[0:-2, 0]) - 
                dt / (2 * rho * dx) * 
                (p[1:-1, 1] - p[1:-1, -1]) + 
                nu * (dt / dx**2 * 
                (un[1:-1, 1] - 2 * un[1:-1, 0] + un[1:-1, -1]) +
                dt / dy**2 *
                (un[2:, 0] - 2 * un[1:-1, 0] + un[0:-2, 0])) + F * dt)

    # Periodic BC v @ x = 2
    v[1:-1, -1] = (vn[1:-1, -1] - un[1:-1, -1] * dt / dx *
                (vn[1:-1, -1] - vn[1:-1, -2]) - 
                vn[1:-1, -1] * dt / dy *
                (vn[1:-1, -1] - vn[0:-2, -1]) -
                dt / (2 * rho * dy) * 
                (p[2:, -1] - p[0:-2, -1]) +
                nu * (dt / dx**2 *
                (vn[1:-1, 0] - 2 * vn[1:-1, -1] + vn[1:-1, -2]) +
                dt / dy**2 *
                (vn[2:, -1] - 2 * vn[1:-1, -1] + vn[0:-2, -1])))

    # Periodic BC v @ x = 0
    v[1:-1, 0] = (vn[1:-1, 0] - un[1:-1, 0] * dt / dx *
                (vn[1:-1, 0] - vn[1:-1, -1]) -
                vn[1:-1, 0] * dt / dy *
                (vn[1:-1, 0] - vn[0:-2, 0]) -
                dt / (2 * rho * dy) * 
                (p[2:, 0] - p[0:-2, 0]) +
                nu * (dt / dx**2 * 
                (vn[1:-1, 1] - 2 * vn[1:-1, 0] + vn[1:-1, -1]) +
                dt / dy**2 * 
                (vn[2:, 0] - 2 * vn[1:-1, 0] + vn[0:-2, 0])))


    # Wall BC: u,v = 0 @ y = 0,2
    u[0, :] = 0
    u[-1, :] = 0
    v[0, :] = 0
    v[-1, :]=0
    
    # udiff = (np.sum(u) - np.sum(un)) / np.sum(u)

    print(stepcount)

    if stepcount % displayfreq == 0:
        # fig = plt.figure(figsize = (11,7), dpi=100)
        # plt.quiver(X[::3, ::3], Y[::3, ::3], u[::3, ::3], v[::3, ::3])

        fig = plt.figure(figsize = (11,7), dpi=100)
        plt.quiver(X, Y, u, v)

        plt.show()

    stepcount += 1