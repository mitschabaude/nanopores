import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, acos, sin
#from aHem_array import *

#x = np.load('exit_x.npy')
#y = np.load('exit_y.npy')
#z = np.load('exit_z.npy')
time = np.load('timer.npy')
counter = np.load('counter.npy')
print float(counter[0])/np.sum(counter)

maxtime=5e6
fineness = 1e4
timestep=maxtime/float(fineness)
exittime = np.linspace(0.,1.,fineness)*maxtime
exitprop = np.zeros(fineness)
for index in range(int(fineness)):
    exitprop[index]=np.where(time<timestep*(index+1))[0].shape[0]
exitprop*=1./np.sum(counter)
plt.semilogx(exittime,exitprop)
plt.scatter(exittime,exitprop)
plt.show()

#def radius(x,y):
#    return sqrt(x**2+y**2)
#def det(a,b,c,d):
#	return a*d-b*c
#
#def normal(ax,ay,bx,by,px,py):
#	AP2=(ax-px)**2+(ay-py)**2
#	BP2=(bx-px)**2+(by-py)**2
#	AB2=(ax-bx)**2+(ay-by)**2
#	AB=sqrt(AB2)
#	c = (AP2-BP2+AB2)/(2*AB)
#	if c>0. and c<AB:
#		if AP2<=c**2:
#			return 0.
#		else:
#			return sqrt(AP2-c**2)
#	else:
#		return 100.
#
#
#def distance_to_surface(rad,z):
#	if z>3.:
#		return 100.
#	elif rad>8.:
#		return 100.
#	else:
#		size=X_aHem.shape[0]
#		D = np.zeros(size)
#		E = np.zeros(size)
#		for index in range(size):
#			D[index] = radius(rad-X_aHem[index][0],z-X_aHem[index][2])
#			E[index] = normal(X_aHem[(index-1)][0],X_aHem[(index-1)][2],X_aHem[(index)][0],X_aHem[(index)][2],rad,z)
#		return min([min(D),min(E)])
#
#
#l=x.shape[0]
#r=np.zeros(l)
#for index in range(l):
#    r[index]=radius(x[index],y[index])
#plt.scatter(r,z,color='red')
#
#
#leftend=max(np.max(r),10.)
#x_mem=np.linspace(X_aHem[18][0],leftend,100)
#y_mem=np.zeros(x_mem.shape[0])+X_aHem[18][2]
#size=X_aHem.shape[0]
#X=np.zeros(size+1)
#Y=np.zeros(size+1)
#for index in range(size):
#	X[index]=X_aHem[index][0]
#	Y[index]=X_aHem[index][2]
#X[size]=X[0]
#Y[size]=Y[0]
#plt.plot(X,Y,linewidth='2',color='blue')
#plt.scatter(X,Y,50,color='blue')
#plt.plot(x_mem,y_mem,color='black',linewidth=1)
#plt.show()
