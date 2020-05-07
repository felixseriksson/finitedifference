import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from datetime import datetime
from os import mkdir

from grid import generategrid
from visualisation import makebox, visualisevelocityvectors, visualisefixedsquare, visualisepressure, savesnap, initdirectory
from boundaryconditions import build_up_b, pressure_poisson_periodic, noslipatwall, noslipatfixedsquare



# variable declarations
nx = 101
ny = 101
nt = 50000
nit = 50 
c = 1
displayfreq = 10

# make grid
gridsizex = 16
gridsizey = 4

dx, dy, X, Y = generategrid(gridsizex, gridsizey, nx, ny)

##physical variables
rho = 1
nu = 0.5
F = 0
dt = 0.001

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

    b = build_up_b(rho, dt, dx, dy, u, v)
    p = pressure_poisson_periodic(p, dx, dy, nit, b, shape=None)

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
    
    # Inlet Condition
    u[1:-1, 0] = 1

    # # Periodic BC u @ x = right
    # u[1:-1, -1] = (un[1:-1, -1] - un[1:-1, -1] * dt / dx * 
    #             (un[1:-1, -1] - un[1:-1, -2]) -
    #             vn[1:-1, -1] * dt / dy * 
    #             (un[1:-1, -1] - un[0:-2, -1]) -
    #             dt / (2 * rho * dx) *
    #             (p[1:-1, 0] - p[1:-1, -2]) + 
    #             nu * (dt / dx**2 * 
    #             (un[1:-1, 0] - 2 * un[1:-1,-1] + un[1:-1, -2]) +
    #             dt / dy**2 * 
    #             (un[2:, -1] - 2 * un[1:-1, -1] + un[0:-2, -1])) + F * dt)

    # # Periodic BC u @ x = 0
    # u[1:-1, 0] = (un[1:-1, 0] - un[1:-1, 0] * dt / dx *
    #             (un[1:-1, 0] - un[1:-1, -1]) -
    #             vn[1:-1, 0] * dt / dy * 
    #             (un[1:-1, 0] - un[0:-2, 0]) - 
    #             dt / (2 * rho * dx) * 
    #             (p[1:-1, 1] - p[1:-1, -1]) + 
    #             nu * (dt / dx**2 * 
    #             (un[1:-1, 1] - 2 * un[1:-1, 0] + un[1:-1, -1]) +
    #             dt / dy**2 *
    #             (un[2:, 0] - 2 * un[1:-1, 0] + un[0:-2, 0])) + F * dt)

    # # Periodic BC v @ x = right
    # v[1:-1, -1] = (vn[1:-1, -1] - un[1:-1, -1] * dt / dx *
    #             (vn[1:-1, -1] - vn[1:-1, -2]) - 
    #             vn[1:-1, -1] * dt / dy *
    #             (vn[1:-1, -1] - vn[0:-2, -1]) -
    #             dt / (2 * rho * dy) * 
    #             (p[2:, -1] - p[0:-2, -1]) +
    #             nu * (dt / dx**2 *
    #             (vn[1:-1, 0] - 2 * vn[1:-1, -1] + vn[1:-1, -2]) +
    #             dt / dy**2 *
    #             (vn[2:, -1] - 2 * vn[1:-1, -1] + vn[0:-2, -1])))

    # # Periodic BC v @ x = 0
    # v[1:-1, 0] = (vn[1:-1, 0] - un[1:-1, 0] * dt / dx *
    #             (vn[1:-1, 0] - vn[1:-1, -1]) -
    #             vn[1:-1, 0] * dt / dy *
    #             (vn[1:-1, 0] - vn[0:-2, 0]) -
    #             dt / (2 * rho * dy) * 
    #             (p[2:, 0] - p[0:-2, 0]) +
    #             nu * (dt / dx**2 * 
    #             (vn[1:-1, 1] - 2 * vn[1:-1, 0] + vn[1:-1, -1]) +
    #             dt / dy**2 * 
    #             (vn[2:, 0] - 2 * vn[1:-1, 0] + vn[0:-2, 0])))
    # v[1:-1, 0] = 0

    # Wall BC: u,v = 0 @ y = bottom, top
    u, v = noslipatwall(u, v)
    u, v = noslipatfixedsquare(u, v)
    # u[25:75, 13:26] = 0
    # v[25:75, 13:26] = 0

    # Obstacle BC: u,v = 0 @ obstacle borders

    print(stepcount)
    if stepcount == 0:
        now = datetime.now()
        dtstring = now.strftime("%d%m%y-%H%M%S")
        dirstring = "C:/Users/felix/Documents/GitHub/gymnasiearbete/model/images/" + str(dtstring)
        mkdir(dirstring)

    if stepcount % displayfreq == 0:

        fig = plt.figure(figsize = (11,7), dpi=100)
        plt.subplot(211)
        makebox(gridsizex, gridsizey)
        visualisefixedsquare()
        visualisevelocityvectors(X, Y, u, v, declutter=True)
        plt.axis("scaled")
        plt.subplot(212)
        makebox(gridsizex, gridsizey)
        visualisefixedsquare()
        visualisepressure(X, Y, p)
        plt.axis("scaled")

        savesnap(dtstring, stepcount)
        plt.close()
        #plt.show()

    stepcount += 1