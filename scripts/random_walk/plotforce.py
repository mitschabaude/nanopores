from math import sqrt
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import numpy as np
from aHem_array_2d import *
from calculateforce import *
F=calculateforce(clscale=10., tol=1e-1)
#plt.ion()
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

leftend=10.
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
axes.set_ylim([-6,5])

bigpath=np.array([[0,0],[0,-5.41],[3,-5.41],[3,0]])
bbPath=mplPath.Path(X_aHem_2d)
bbPath_2=mplPath.Path(bigpath)
def plot_arrows(fac=0.1,p1=[-10.,10.],p2=[5.,-5.],Nx=50,Ny=40,bbPath=bbPath):
    x_vec=np.linspace(p1[0],p1[1],Nx)
    z_vec=np.linspace(p2[0],p2[1],Ny)
    plt.plot(X,Y,linewidth='2',color='blue')
    plt.plot(X_2,Y,linewidth=2,color='blue')
    plt.plot(x_mem,y_mem,color='black',linewidth=1)
    plt.plot(x_mem_2,y_mem,color='black',linewidth=1)
    for x in x_vec:
        for z in z_vec:
            if not bbPath.contains_point((x,z)) and not bbPath.contains_point((-x,z)) and not (z<0. and x>-3. and x<3.):
                Force=F(argument(x,0,z))
                Fx=Force[0]*1e12
                Fz=Force[2]*1e12
                ax.arrow(x,z,fac*Fx,fac*Fz/2.,head_width=fac*0.3*abs(Fz), head_length=fac*abs(Fz/2.),linewidth=2)
plot_arrows() # only mouth of nanopore
plt.show()
