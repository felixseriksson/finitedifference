
import numpy as np
from matplotlib import pyplot as plt
import time
import sys

nx = 41
dx = 2/(nx -1)
nt = 100
dt = 0.025
c = 1

u = np.ones(nx)
u[int(0.5/dx):int(1/dx+1)] = 2
print(u)

plt.plot(np.linspace(0, 2, nx), u)
plt.show()

un = np.ones(nx)

for n in range(nt):
    un = u.copy()
    for i in range(1, nx):
        u[i] = un[i] - c* (dt/dx)* (un[i]- un[i-1])
    if n % 10 == 0:
        plt.plot(np.linspace(0, 2, nx), u)
        plt.show()

plt.plot(np.linspace(0, 2, nx), u)
plt.show()
'''
#lånar denna för en rekursiv sak i Matematik 1
def G(n):
    if n == 0:
        return 1
    elif n == 1:
        return 1
    else:
        return 3*G(n-1) + 3*G(n-2)

print(G(30))
'''