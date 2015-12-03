from math import pi, sqrt
# from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pylab
import sys
from aHem_array import *

counter=np.load('counter.npy') # ahem, molecule, poretop, membrane, bulk
counter[1]=0
hits=np.sum(counter)


exit_rad=np.load('exit_rad.npy')
exit_z=np.load('exit_z_2d.npy')

hist=np.histogram(exit_rad,bins=30,range=(0,40))
x=hist[1]
y=hist[0]
x=np.delete(x,0)
plt.plot(x,y)
plt.show()	