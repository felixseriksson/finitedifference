import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime
from os import mkdir

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

def visualisepressure(X, Y, p):
    plt.pcolormesh(X, Y, p, cmap="autumn")

def visualiseobstacle(shape):
    return 0

def visualisefixedsquare():
    plt.plot([2.1, 3.5],[2.7, 2.7],"k")
    plt.plot([2.1, 3.5],[1.3, 1.3],"k")
    plt.plot([2.1, 2.1],[1.3, 2.7],"k")
    plt.plot([3.5, 3.5],[1.3, 2.7],"k")

def savesnap(dtstring,number):
    plt.savefig("C:\\Users\\felix\\Documents\\GitHub\\gymnasiearbete\\images\\"+str(dtstring)+"\\"+str(number)+".png")

def initdirectory():
    now = datetime.now()
    dtstring = now.strftime("%d%m%y-%H%M%S")
    dirstring = "C:\\Users\\felix\\Documents\\GitHub\\gymnasiearbete\\images\\" + str(dtstring)
    mkdir(dirstring)