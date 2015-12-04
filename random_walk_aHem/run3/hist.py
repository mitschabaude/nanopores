from math import pi, sqrt
# from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pylab as P
import sys
from aHem_array import *

counter=np.load('counter.npy') # ahem, molecule, poretop, membrane, bulk
counter[1]=0
hits=float(np.sum(counter))
good=float(counter[0]+counter[2]+counter[3])

bin1=np.linspace(0,200,50)
bin2=np.linspace(0,20,50)
bin3=np.linspace(0,5,20)
bin=np.concatenate((bin1,bin2))
exit_rad=np.load('exit_rad.npy')
exit_z=np.load('exit_z_2d.npy')

exit=0.

for index in range(exit_rad.shape[0]):
	if exit_rad[index]<=2.16:
		exit+=1.
print exit/hits


# hist=np.histogram(exit_rad,bins=bin)
# x=hist[1]
# y=hist[0]
# x=np.delete(x,0)
# plt.scatter(x,y)
# plt.show()	

n, bins, patches = P.hist(exit_rad, bin2,normed=1, histtype='stepfilled')
P.setp(patches, 'facecolor', 'g', 'alpha', .75)
P.show()