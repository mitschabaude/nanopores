import random,math
from math import sqrt, pi
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pylab

from params_geo import *


mpl.rcParams['legend.fontsize']=10

fig = plt.figure()
#ax = fig.gca(projection='3d')
ax = fig.add_subplot(111, projection='3d')
#plt.hold(True)

X=np.load('X.npy')
Y=np.load('Y.npy')
Z=np.load('Z.npy')
#X_C=np.load('X_C.npy')
#Y_C=np.load('Y_C.npy')
#Z_C=np.load('Z_C.npy')
X=X[:np.floor(X.shape[0]/2.0)-0]
Y=Y[:np.floor(Y.shape[0]/2.0)-0]
Z=Z[:np.floor(Z.shape[0]/2.0)-0]

#path
ax.plot(X,Y,Z,c='green')

#Collisions
#ax.scatter(X_C,Y_C,Z_C,c='red')

#################
show_all=False
#innercylinder
alphb=0.0
#rest
alpha=0.0

# innerCylinder

x=np.linspace(-r0*1e9, r0*1e9, 100)
z=np.linspace(-l0/2.0*1e9, l0/2.0*1e9, 100)
Xic, Zic=np.meshgrid(x, z)
Yic = np.sqrt((r0*1e9)**2-Xic**2)

# Draw parameters
rstride = 40
cstride = 40
ax.plot_surface(Xic, Yic, Zic, alpha=alphb, rstride=rstride, cstride=cstride)
if show_all:
    ax.plot_surface(Xic, -Yic, Zic, alpha=alphb, rstride=rstride, cstride=cstride)

ax.set_xlabel("X")


# outerupperCylinder

x=np.linspace(-r1*1e9, r1*1e9, 100)
z=np.linspace(l1/2.0*1e9, l0/2.0*1e9, 100)
Xouc, Zouc=np.meshgrid(x, z)
Youc = np.sqrt((r1*1e9)**2-Xouc**2)

# Draw parameters
rstride = 50
cstride = 50
ax.plot_surface(Xouc, Youc, Zouc, alpha=alpha, rstride=rstride, cstride=cstride)
if show_all:
    ax.plot_surface(Xouc, -Youc, Zouc, alpha=alpha, rstride=rstride, cstride=cstride)



# outerlowerCylinder

x=np.linspace(-r1*1e9, r1*1e9, 100)
z=np.linspace(-l0/2.0*1e9, -l1/2.0*1e9, 100)
Xodc, Zodc=np.meshgrid(x, z)
Yodc = np.sqrt((r1*1e9)**2-Xodc**2)

# Draw parameters
#rstride = 40
#cstride = 40
ax.plot_surface(Xodc, Yodc, Zodc, alpha=alpha, rstride=rstride, cstride=cstride)
if show_all:
    ax.plot_surface(Xodc, -Yodc, Zodc, alpha=alpha, rstride=rstride, cstride=cstride)


# upperTop

u = np.linspace(0, np.pi, 100)
v = np.linspace(r0*1e9, r1*1e9, 10)

x = np.outer(v, np.cos(u))
y = np.outer(v, np.sin(u))
z = np.zeros(shape=x.shape)+l0/2.0*1e9
ax.plot_surface(x, y, z,  rstride=15, cstride=15, color='b', alpha=alpha)
if show_all:
    ax.plot_surface(x, -y, z,  rstride=15, cstride=15, color='b', alpha=alpha)

# lowerTop

u = np.linspace(0, np.pi, 100)
v = np.linspace(r0*1e9, r1*1e9, 10)

x = np.outer(v, np.cos(u))
y = np.outer(v, np.sin(u))
z = np.zeros(shape=x.shape)-l0/2.0*1e9
ax.plot_surface(x, y, z,  rstride=15, cstride=15, color='b', alpha=alpha)
if show_all:
    ax.plot_surface(x, -y, z,  rstride=15, cstride=15, color='b', alpha=alpha)

# upperMembran

u = np.linspace(0, np.pi, 100)
v = np.linspace(r1*1e9, Rz*1e9, 10)

x = np.outer(v, np.cos(u))
y = np.outer(v, np.sin(u))
z = np.zeros(shape=x.shape)+l1/2.0*1e9
ax.plot_surface(x, y, z,  rstride=15, cstride=15, color='b', alpha=alpha)
if show_all:
    ax.plot_surface(x, -y, z,  rstride=15, cstride=15, color='b', alpha=alpha)

# lowerMembran

u = np.linspace(0, np.pi, 100)
v = np.linspace(r1*1e9, Rz*1e9, 10)

x = np.outer(v, np.cos(u))
y = np.outer(v, np.sin(u))
z = np.zeros(shape=x.shape)-l1/2.0*1e9
ax.plot_surface(x, y, z,  rstride=15, cstride=15, color='b', alpha=alpha)
if show_all:
    ax.plot_surface(x, -y, z,  rstride=15, cstride=15, color='b', alpha=alpha)
#################

# Sphere
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = rMolecule*1e9*np.outer(np.cos(u), np.sin(v))
y = rMolecule*1e9*np.outer(np.sin(u), np.sin(v))
z = rMolecule*1e9*np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x+X[X.shape[0]-1], y+Y[Y.shape[0]-1], z+Z[Z.shape[0]-1],  rstride=12, cstride=5, color='r',alpha=1.0)

ax.scatter([-10,10,-10,10,-10,10,-10,10],[10,10,-10,-10,10,10,-10,-10],[-10,-10,-10,-10,10,10,10,10],c='w',alpha=0.0)

plt.show()
