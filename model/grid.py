import numpy as np

def generategrid(xsize, ysize, xsteps, ysteps):
    dx = xsize / (xsteps - 1)
    dy = ysize / (ysteps - 1)
    
    x = np.linspace(0, xsize, xsteps)
    y = np.linspace(0, ysize, ysteps)

    meshX, meshY = np.meshgrid(x, y)
    return dx, dy, meshX, meshY