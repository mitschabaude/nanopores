from math import pi
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pylab
import sys

fig = plt.figure()
ax=Axes3D(fig)

EXIT_X=np.load('exit_x_1.npy')
EXIT_Y=np.load('exit_y_1.npy')
EXIT_Z=np.load('exit_z_1.npy')
counter=np.load('counter_1.npy')
hits = float(np.sum(counter))
print
print
print 'PROBABILITY REACHING:'
print 
print 'PAC-MAN: ', float(counter[1])/hits
print
print 'MEMBRANE: ', float(counter[3])/hits
print
print 'A-HEMOLYSIN: ', float(counter[0])/hits
print
print 'EXIT: ', float(counter[2])/hits
print
print 'INFINITY: ', float(counter[4])/hits
print

delete=[]
for index in range(EXIT_X.shape[0]):
    if EXIT_Z[index]<-5:
        delete.append(index)

for index in reversed(delete):
    EXIT_X = np.delete(EXIT_X, index)
    EXIT_Y = np.delete(EXIT_Y, index)
    EXIT_Z = np.delete(EXIT_Z, index)

delete=[]
for index in range(EXIT_X.shape[0]):
    if ((EXIT_X[index]-5)**2+(EXIT_Y[index])**2+(EXIT_Z[index]-10)**2)<2.2**2:
        delete.append(index)

for index in reversed(delete):
    EXIT_X = np.delete(EXIT_X, index)
    EXIT_Y = np.delete(EXIT_Y, index)
    EXIT_Z = np.delete(EXIT_Z, index)

# delete=[]
# for index in range(EXIT_X.shape[0]):
#     if EXIT_Y[index]<0:
#         delete.append(index)

# for index in reversed(delete):
#     EXIT_X = np.delete(EXIT_X, index)
#     EXIT_Y = np.delete(EXIT_Y, index)
#     EXIT_Z = np.delete(EXIT_Z, index)


show_half=False
fac=2
if show_half:
    fac=1
from aHem_array import *

x_aHem, z_aHem = np.zeros_like(X_aHem), np.zeros_like(X_aHem)
for index in range(x_aHem.shape[0]):
    x_aHem[index], z_aHem[index] = X_aHem[index][0], X_aHem[index][2]
x_aHem = np.append(x_aHem,x_aHem[0])
z_aHem = np.append(z_aHem,z_aHem[0])
angle = np.linspace(0,fac*pi,20)
theta, mesh_x = np.meshgrid(angle,x_aHem)
theta, mesh_z = np.meshgrid(angle,z_aHem)

X_plot = mesh_x*np.cos(theta)
Y_plot = mesh_x*np.sin(theta)
Z_plot = mesh_z
ax.plot_wireframe(X_plot, Y_plot, Z_plot,color='blue', rstride = 1, cstride = 1, alpha=.1)

if False:#show_half:
    x=x_aHem.tolist()
    left=-x_aHem
    x_=left.tolist()
    y=np.zeros_like(x_aHem).tolist()
    z=z_aHem.tolist()
    verts1=[zip(x,y,z)]
    verts2=[zip(x_,y,z)]
    ax.add_collection3d(Poly3DCollection(verts1,color='blue',alpha=1.0))
    ax.add_collection3d(Poly3DCollection(verts2,color='blue',alpha=1.0))
ax.set_xlim3d(-10,10)
ax.set_ylim3d(-10,10)
ax.set_zlim3d(-10,10)
# Sphere
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

# for index in range(EXIT_X.shape[0]-1):
#     x = 0.1*np.outer(np.cos(u), np.sin(v))
#     y = 0.1*np.outer(np.sin(u), np.sin(v))
#     z = 0.1*np.outer(np.ones(np.size(u)), np.cos(v))
#     ax.plot_surface(x+EXIT_X[index], y+EXIT_Y[index], z+EXIT_Z[index],  rstride=12, cstride=5, color='r',alpha=1.0)
ax.scatter(EXIT_X,EXIT_Y,EXIT_Z,color='red')


x1 = 2.2*np.outer(np.cos(u), np.sin(v))
y1 = 2.2*np.outer(np.sin(u), np.sin(v))
z1 = 2.2*np.outer(np.ones(np.size(u)), np.cos(v))
# ax.plot_wireframe(x1+5., y1+0., z1+10., rstride=12, cstride=5, color='blue',alpha=1.0)
plt.show()
