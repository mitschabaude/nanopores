import colormaps
from math import sqrt
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import numpy as np
from aHem_array_2d import *
from calculateforce import loadforces
F, Fel, Fdrag = loadforces()
#plt.ion()
fig=plt.figure(figsize=(18,12))
ax=plt.axes()
def argument(x,y,z):
    return np.array([float(x),float(y),float(z)])
#def F(vec):
#    return [0.,0.,-2e-13]
def radius(x,y):
    return sqrt(x**2+y**2)
def sgn(x):
    if x<0:
        return -1
    elif x>0:
        return 1
    else:
        return 0

leftend=15.
x_mem=np.linspace(X_aHem_2d[18][0],leftend,100)
y_mem=np.zeros(x_mem.shape[0])+X_aHem_2d[18][1]
x_mem_2=-x_mem
size=X_aHem_2d.shape[0]
X=np.zeros(size+1)
Y=np.zeros(size+1)
for index in range(size):
	X[index]=X_aHem_2d[index][0]
	Y[index]=X_aHem_2d[index][1]
X[size]=X[0]
Y[size]=Y[0]
X_2=-X

# whole domain: fac=0.1,p2=[

axes=plt.gca()
axes.set_ylim([-10,15])
axes.set_xlim([-15,15])

bbPath=mplPath.Path(X_aHem_2d)
Ny = 40
Nx = 50
plt.plot(X,Y,linewidth='2',color='black')
plt.plot(X_2,Y,linewidth=2,color='black')
plt.plot(x_mem,y_mem,color='black',linewidth=1)
plt.plot(x_mem_2,y_mem,color='black',linewidth=1)

Y, X = np.mgrid[-5:15:Ny*1j, -15:15:Nx*1j]
U = np.zeros((Ny,Nx))
V = np.zeros((Ny,Nx))
for y in range(Ny):
    for x in range(Nx):
        if bbPath.contains_point((X[y][x],Y[y][x])) or bbPath.contains_point((-X[y][x],Y[y][x])):
            U[y][x] = 0
            V[y][x] = 0
        else:
            F=Fel(argument(X[y][x],0,Y[y][x]))
            U[y][x] = F[0]
            V[y][x] = F[2]

plt.streamplot(X,Y,U,V,arrowsize=2, linewidth=1., density=4., cmap=colormaps.viridis, color=np.log10(np.sqrt(U*U+V*V)))
#plt.savefig('el.eps')
plt.show()
