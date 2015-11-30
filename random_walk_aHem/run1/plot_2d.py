from math import pi, sqrt
# from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pylab
import sys
from aHem_array import *

exit_rad=np.load('exit_rad.npy')
exit_z=np.load('exit_z_2d.npy')

x_aHem, z_aHem = np.zeros_like(X_aHem), np.zeros_like(X_aHem)
for index in range(x_aHem.shape[0]):
    x_aHem[index], z_aHem[index] = X_aHem[index][0], X_aHem[index][2]
x_aHem = np.append(x_aHem,x_aHem[0])
z_aHem = np.append(z_aHem,z_aHem[0])

x=np.linspace(4.96,15)
y=np.zeros_like(x)-5.41

plt.plot(x,y)
plt.plot(exit_rad,exit_z,'o')
plt.axis([0, 15, -10, 1])
plt.plot(x_aHem,z_aHem)
plt.show()