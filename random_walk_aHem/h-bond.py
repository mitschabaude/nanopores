from math import sqrt, pi
import math
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from math import exp
D_kcal = 100 # kcal/mol
n = 55 # 1/nm
r_e = 1. #nm - equilibrium distance
D_en = D_kcal*4184e9 # from kcal to N*nm
mol = 6.022e-23 # 1 atom

h=0.01
a=0.5
b=2.
a1=.1
b1=5.5
X=np.arange(a,b,h)
X1=np.arange(a1,b1,h)
size1=X1.shape[0]
size=X.shape[0]

def V(r,r_e,n,D,mol):
    if r>0.2 and r<3.:
        d_r = r-r_e
        return D*mol*(1-math.exp(-n*(d_r**2)/(2*r)))
    else:
        return D*mol
def Vdiff(r,r_e,n,D,mol):
    if r>0.2 and r<3.:
        d_r = r-r_e
        return -D_en*mol*exp(-n*(d_r**2)/(2*r))*(n*d_r**2/(2*r**2)-n*d_r/r)
    else:
        return 0.
def surf(r):
    return 1e-10/(r**12)
def surfdiff(r):
    return -12*1e-10/(r**13)

C=np.zeros(size)

Y=np.zeros(size)
for index in range(size):
	Y[index]=V(X[index],r_e,n,D_en,mol)
	C[index]=surf(X[index])

P=C+Y


D_kcal = 100 # kcal/mol
n = 55 # 1/nm
r_e = 1. #nm - equilibrium distance
D_en = D_kcal*4184e9 # from kcal to N*nm
mol = 6.022e-23 # 1 atom

F = np.zeros(Y.shape[0])
F1 = np.zeros(size)
for index in range(1,size-1):
    F[index]=-(P[index+1]-P[index-1])/(2*h)
    d = (b-a)/size*(index)+a
    F1[index] = -Vdiff(d,r_e,n,D_en,mol)-surfdiff(d)

F[0]=F[1]
F[size-1]=F[size-2]

red_patch = mpatches.Patch(color='red', label='Sum of potentials')
green_patch = mpatches.Patch(color='green', label='H-Bond potential')
blue_patch = mpatches.Patch(color='blue', label='Surface potential')
orange_patch = mpatches.Patch(color='orange', label='-Vdiff-surfdiff')
yellow_patch = mpatches.Patch(color='yellow', label='-(P[+1]-P[-1]/2h')

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
plt.ylim([-1e-7,1.5e-7])
plt.legend(handles=[orange_patch,yellow_patch])
plt.plot(X,F,color='yellow',linewidth='2')
plt.plot(X,np.zeros(size),color='black',linewidth='2')
plt.plot(X,F1,color='orange',linewidth='2')
plt.show()
