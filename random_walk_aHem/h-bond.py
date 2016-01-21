from math import sqrt, pi
import math
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

D_kcal = 100 # kcal/mol
n = 55 # 1/nm
r_e = 1. #nm - equilibrium distance
D = D_kcal*4184e9 # from kcal to N*nm
mol = 6.022e-23 # 1 atom

h=0.01
X=np.arange(0.05,3.,h)
size=X.shape[0]

def V(r,r_e,n,D,mol):
	d_r = r-r_e
	return D*mol*(1-math.exp(-n*(d_r**2)/(2*r)))

def logi(r,a):
	return math.exp(-r*a)/(1+math.exp(-r*a))
#	return 1.
def lin(r,a):
	if r<-1./a:
		return 1
	elif r>1./a:
		return 0
	else:
		return 0.5-2./a*r


def surf(r):
	return 1e-5/((r+1.4)**12)

C=np.zeros(size)

Y=np.zeros(size)
for index in range(size):
	Y[index]=V(X[index],r_e,n,D,mol)
	C[index]=surf(X[index])

P=C+Y

F = np.zeros(Y.shape[0])
for index in range(1,size-1):
	F[index]=-(P[index+1]-P[index-1])/(2*h)
F[0]=F[1]
F[size-1]=F[size-2]
print 'Sum of Potential = 1e-5*(1/(x+1.4))^12'
red_patch = mpatches.Patch(color='red', label='Sum of potentials')
green_patch = mpatches.Patch(color='green', label='H-Bond potential')
blue_patch = mpatches.Patch(color='blue', label='Surface potential')
orange_patch = mpatches.Patch(color='orange', label='Force')

ax = plt.gca().add_artist(plt.legend(handles=[blue_patch],loc=1))
plt.legend(handles=[green_patch],loc=4)
plt.plot(X,Y,color='green',linewidth='2')
plt.plot(X,C,color='blue',linewidth='2')
plt.ylim([-0.1e-7,1.2e-7])
plt.show()
ax = plt.gca().add_artist(plt.legend(handles=[blue_patch],loc=1))
ax = plt.gca().add_artist(plt.legend(handles=[green_patch],loc=4))
plt.legend(handles=[red_patch],loc=3)
plt.plot(X,Y,color='green',linewidth='1')
plt.plot(X,C,color='blue',linewidth='1')
plt.plot(X,P,color='red',linewidth='2')
plt.ylim([-0.1e-7,1.2e-7])
plt.show()
plt.legend(handles=[orange_patch])
plt.plot(X,F,color='red',linewidth='2')
plt.plot(X,np.zeros(size),color='black',linewidth='2')
plt.show()
