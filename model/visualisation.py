from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def makebox(gridsizex, gridsizey):
    plt.plot([0,gridsizex],[gridsizey, gridsizey],"k")
    plt.plot([0,gridsizex],[0, 0],"k")
    plt.plot([0,0],[0, gridsizey],"k--")
    plt.plot([gridsizex, gridsizex],[0, gridsizey],"k:")


def visualisevelocityvectors(X, Y, u, v, declutter=False):
    if declutter:
        plt.quiver(X[::3, ::3], Y[::3, ::3], u[::3, ::3], v[::3, ::3])
    else:
        plt.quiver(X, Y, u, v)