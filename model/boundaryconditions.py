import numpy as np

from obstacles import *

def build_up_b(rho, dt, dx, dy, u, v):
    b = np.zeros_like(u)
    b[1:-1, 1:-1] = (rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) +
                                      (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
                            ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 -
                            2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
                                 (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx))-
                            ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2))
    
    # Periodic BC Pressure @ x = right
    b[1:-1, -1] = (rho * (1 / dt * ((u[1:-1, 0] - u[1:-1,-2]) / (2 * dx) +
                                    (v[2:, -1] - v[0:-2, -1]) / (2 * dy)) -
                          ((u[1:-1, 0] - u[1:-1, -2]) / (2 * dx))**2 -
                          2 * ((u[2:, -1] - u[0:-2, -1]) / (2 * dy) *
                               (v[1:-1, 0] - v[1:-1, -2]) / (2 * dx)) -
                          ((v[2:, -1] - v[0:-2, -1]) / (2 * dy))**2))

    # Periodic BC Pressure @ x = 0
    b[1:-1, 0] = (rho * (1 / dt * ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx) +
                                   (v[2:, 0] - v[0:-2, 0]) / (2 * dy)) -
                         ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx))**2 -
                         2 * ((u[2:, 0] - u[0:-2, 0]) / (2 * dy) *
                              (v[1:-1, 1] - v[1:-1, -1]) / (2 * dx))-
                         ((v[2:, 0] - v[0:-2, 0]) / (2 * dy))**2))
    
    return b

def pressure_poisson_periodic(p, dx, dy, nit, b, shape=None):
    pn = np.empty_like(p)
    
    for _ in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 +
                          (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) /
                         (2 * (dx**2 + dy**2)) -
                         dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, 1:-1])

        # Periodic BC Pressure @ x = right
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
        p = pressureatwall(p)

        # pressure boundary conditions for obstacle
        # obstacledict[shape]()
        p = pressureatfixedsquare(p)
    
    return p

def noslipatwall(u,v):
    u[0, :] = 0
    u[-1, :] = 0
    v[0, :] = 0
    v[-1, :] = 0
    return u, v

def pressureatwall(p):
    p[-1, :] = p[-2, :]  # dp/dy = 0 at y = top
    p[0, :] = p[1, :]  # dp/dy = 0 at y = 0
    return p


# no-slip function for obstacles
def noslipatobstacle(shape):
    return 0

def noslipatsquare(u, v, sidelength=1.4, centerx=2.7, centery=2):
    # default: hörn i (2, 1.3), (3.4, 1.3), (2, 2.7) och (3.4, 2.7)
    # (anpassat för medium i 4 x 16 ruta)
    return 0

def noslipatfixedsquare(u, v):
    # this function has no variable input and works ONLY for
    # a 4 x 16 box with nx = ny = 101
    for indexy in range(32, 69):
        for indexx in range(12, 23):
            u[indexy, indexx] = 0
            v[indexy, indexx] = 0
    return u, v

# pressure bc function for obstacles, called in pressure_poisson_periodic
def pressureatobstacle(shape):
    return 0

def pressureatfixedsquare(p):
    # this function has no variable input and works ONLY for
    # a 4 x 16 box with nx = ny = 101
    for indexy in range(32, 69):
        for indexx in range(12, 23):
            p[indexy, indexx] = 0
    for indexx in range(12, 23):
        p[32, indexx] = p[31, indexx]
        p[68, indexx] = p[69, indexx]
    for indexy in range(32, 69):
        p[indexy, 12] = p[indexy, 11]
        p[indexy, 22] = p[indexy, 23]
    return p