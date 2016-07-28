from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
#import colormaps
from matplotlib import cm
from math import sqrt
import matplotlib.pyplot as plt
import matplotlib.path as mplPath
import numpy as np
from calculateforce import loadforces
Fel_, Fdrag_ = loadforces()
X_How_2d = np.array([[1.,4.5],[2.5,4.5],[2.5,1.1],[10.,1.1],[10.,-1.1],[2.5,-1.1],[2.5,-4.5],[1.,-4.5]])
def Fel(x,y,z):
    if x==0. and y==0.:
        return [0.,0.,Fel_(np.array([0,z]))[1]]
    else:
        rad=sqrt(x**2+y**2)
        x0=x/rad
        y0=y/rad
        return [x0*Fel_(rad,z)[0],y0*Fel_(rad,z)[0],Fel_(rad,z)[1]]
def Fdrag(x,y,z):
    if x==0. and y==0.:
        return [0.,0.,Fdrag_(np.array([0,z]))[1]]
    else:
        rad=sqrt(x**2+y**2)
        x0=x/rad
        y0=y/rad
        return [x0*Fdrag_(rad,z)[0],y0*Fdrag_(rad,z)[0],Fdrag_(rad,z)[1]]

#plt.ion()
#fig1=plt.figure(figsize=(18,12))
#fig=fig1.add_subplot()
#bar=fig1.add_subplot()
bar, fig = plt.subplots(figsize=(12,8))
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
#x_mem=np.linspace(X_How_2d[18][0],leftend,100)
#y_mem=np.zeros(x_mem.shape[0])+X_How_2d[18][1]
#x_mem_2=-x_mem
size=X_How_2d.shape[0]
X=np.zeros(size+1)
Y=np.zeros(size+1)
for index in range(size):
	X[index]=X_How_2d[index][0]
	Y[index]=X_How_2d[index][1]
X[size]=X[0]
Y[size]=Y[0]
X_2=-X

# whole domain: fac=0.1,p2=[

axes=plt.gca()
axes.set_ylim([-5,10])
axes.set_xlim([-10,10])

bbPath=mplPath.Path(X_How_2d)
Ny = 40
Nx = 50
plt.plot(X,Y,linewidth=3,color='black')
plt.plot(X_2,Y,linewidth=3,color='black')
#plt.plot(x_mem,y_mem,color='black',linewidth=1)
#plt.plot(x_mem_2,y_mem,color='black',linewidth=1)

Y, X = np.mgrid[-5:10:Ny*1j, -10:10:Nx*1j] #-5:30 and -30:30
U = np.zeros((Ny,Nx))
V = np.zeros((Ny,Nx))
for y in range(Ny):
    for x in range(Nx):
        if bbPath.contains_point((X[y][x],Y[y][x])) or bbPath.contains_point((-X[y][x],Y[y][x])):
            U[y][x] = 0
            V[y][x] = 0
        else:
            F=Fel(X[y][x],0.,Y[y][x])
            if sqrt(F[0]**2+F[1]**2+F[2]**2)<1e-16:
                U[y][x] = 0.
                V[y][x] = 0.
            else:
                U[y][x] = F[0]
                V[y][x] = F[2]
#for y in range(Ny):
#    for x in range(Nx):
#        if bbPath.contains_point((X[y][x],Y[y][x])) or bbPath.contains_point((-X[y][x],Y[y][x])):
#            U[y][x] = 0
#            V[y][x] = 0
#        else:
#            if Y[y][x]<-5.5 and (X[y][x]<-1.2 or X[y][x]>1.2):
#                U[y][x] = 0
#                V[y][x] = 0
#            else:
#                F=Fdrag(argument(X[y][x],0,Y[y][x]))
#                U[y][x] = F[0]
#                V[y][x] = F[2]

strm = plt.streamplot(X,Y,U,V,arrowsize=3, linewidth=2., density=2.5, cmap=cm.viridis, color=np.log10(np.sqrt(U*U+V*V)))
bar.colorbar(strm.lines)
ax.set_xlabel('X coordinate [nm]', fontsize=20)
ax.set_ylabel('Z coordinate [nm]', fontsize=20)
plt.title(r'$F_{\mathrm{el}}(x)$', fontsize=25)
#plt.show()
plt.savefig('el.eps')
