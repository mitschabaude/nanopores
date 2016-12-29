import sys
from math import sqrt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
#from numpy.random import random
import matplotlib.pyplot as plt
from folders import fields
import nanopores.geometries.pughpore as pughpore
import nanopores
from sobol.sobol_seq import i4_sobol_generate as sobol

up = nanopores.user_params(pughpore.params, k=3)

R = up.R
H = up.H
l0 = up.l0
l1 = up.l1
l2 = up.l2
l3 = up.l3
l4 = up.l4
hpore = up.hpore
hmem = up.hmem
h2 = up.h2
h1 = up.h1
h4 = up.h4
rMolecule = up.rMolecule
eps = 0.1
r = rMolecule + eps
p0=hpore/2.
p1=p0-h1
p2=p0-h2
p3=-hpore/2.

def fac(z):
    if z>=p3 and z<=p2:
        return l3/2.-r
    elif z>p2 and z<p2+r:
        x=z-p2
        return -sqrt(r**2-x**2)+l3/2.
    elif z>=p2+r and z<=p1:
        return l2/2.-r
    elif z>p1 and z<p1+r:
        x=z-p1
        return -sqrt(r**2-x**2)+l2/2.
    elif z>=p1+r and z<=p0:
        return l0/2.-r
    elif z>p0 and z<p0+r:
        x=z-p0
        return -sqrt(r**2-x**2)+l0/2.
    elif z<p3 and z>p3-r:
        x=z-p3
        return -sqrt(r**2-x**2)+l3/2.
    else: return R/2.


k0=3
n=5
#k=up.k
k=0
z=sobol(1,2**n,0)[0]
z=(hpore+30.)*z-hpore/2.-10.
factors=np.array([fac(x) for x in z])

#def R_(z):
#    if z>=p3 and z<=p2:
#        return l3/2.
#    elif z>=p2 and z<=p1:
#        return l2/2.
#    elif z>=p1 and z<=p0:
#        return l0/2.
#    else: return R/2.
#
#X=np.linspace(-H/2.,H/2.,800)
#Y=np.array([fac(X[i]) for i in range(X.shape[0])])
#plt.plot(X,Y)
#plt.plot(X,np.array([R_(X[i]) for i in range(X.shape[0])]),color='r')
#plt.scatter(z,np.zeros(2**n))
#plt.scatter(z,factors)
#plt.show()


XY = sobol(2,2**k0,0)
X_points,Y_points = XY[0],XY[1]
for i in range(k+1): # additional points
    XY = sobol(2,2**(i+k0),2**(i+k0))
    X_points = np.append(X_points,XY[0])
    Y_points = np.append(Y_points,XY[1])

for i in list(reversed(range(X_points.shape[0]))): # cut off other triangle
    if Y_points[i]>X_points[i]:
        X_points = np.delete(X_points,i)
        Y_points = np.delete(Y_points,i)

print '# points in plane = %d\n# z-values =%d\n# totals points= %d'%(
    X_points.shape[0],z.shape[0],X_points.shape[0]*z.shape[0])

X, Y, Z = np.array([]), np.array([]), np.array([]) # concatenate all arrays
for j in range(z.shape[0]):
    Z_p=np.zeros(X_points.shape[0])+z[j]
    X_p = X_points*factors[j]
    Y_p = Y_points*factors[j]
    X = np.append(X,X_p)
    Y = np.append(Y,Y_p)
    Z = np.append(Z,Z_p)
array=[[X[i],Y[i],Z[i]] for i in range(X.shape[0])]

if __name__ == "__main__":
    fields.save_entries("pughx_new", dict(up), x=array, N=len(array))
    fields.update()

    def surfx(y1,y2,z1,z2,d,size,rs,cs):
        Y = np.linspace(y1,y2,size)
        Z = np.linspace(z1,z2,size)
        Y, Z = np.meshgrid(Y,Z)
        X = np.zeros(size)+d
        surf = ax.plot_surface(X,Y,Z, rstride=rs, cstride=cs,
                               alpha=alpha,color=color)

    def surfy(x1,x2,z1,z2,d,size,rs,cs):
        X = np.linspace(x1,x2,size)
        Z = np.linspace(z1,z2,size)
        X, Z = np.meshgrid(X,Z)
        Y = np.zeros(size)+d
        surf = ax.plot_surface(X,Y,Z, rstride=rs, cstride=cs,
                               alpha=alpha,color=color)

    def surfz(x1,x2,y1,y2,d,size,rs,cs):
        X = np.linspace(x1,x2,size)
        Y = np.linspace(y1,y2,size)
        X, Y = np.meshgrid(X,Y)
        Z = np.zeros(size)+d
        surf = ax.plot_surface(X,Y,Z, rstride=rs, cstride=cs,
                               alpha=alpha,color=color)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.view_init(elev=0,azim=270)
    ax.set_xlim([-.5*R,.5*R])
    ax.set_ylim([-.5*R,.5*R])
    ax.set_zlim([-.5*H,.5*H])
    #ax.set_aspect(1)

    size=10
    alpha=.1
    rs, cs = 1, 1
    color='blue'

    #front
    surfy(-.5*l3,.5*l3,-.5*hpore,.5*hpore-h2,.5*l3,size,1,1)
    surfy(-.5*l2,.5*l2,.5*hpore-h2,.5*hpore-h1,.5*l2,size,5,1)
    surfy(-.5*l1,.5*l1,.5*hpore-h1,.5*hpore,.5*l1,size,10,1)
    surfy(.5*l0,-.5*l0,-.5*hpore+hmem,.5*hpore,.5*l0,size,5,5)
    #front-right
    surfy(.5*l3,.5*l2,-.5*hpore,.5*hpore-h2,0.,size,5,5)
    surfy(.5*l2,.5*l1,-.5*hpore,.5*hpore-h1,0.,size,5,5)
    surfy(.5*l1,.5*l0,-.5*hpore+hmem,.5*hpore,0.,size,5,5)
    surfy(.5*l4,.5*R,-.5*hpore,-.5*hpore+hmem,0.,size,5,5)
    #front-left
    surfy(-.5*l3,-.5*l2,-.5*hpore,.5*hpore-h2,0.,size,5,5)
    surfy(-.5*l2,-.5*l1,-.5*hpore,.5*hpore-h1,0.,size,5,5)
    surfy(-.5*l1,-.5*l0,-.5*hpore+hmem,.5*hpore,0.,size,5,5)
    surfy(-.5*l4,-.5*R,-.5*hpore,-.5*hpore+hmem,0.,size,5,5)

    #top-front
    surfz(-.5*l0,.5*l0,.5*l1,.5*l0,.5*hpore,size,10,1)
    surfz(-.5*l1,.5*l1,.5*l2,.5*l1,.5*hpore-h1,size,10,1)
    surfz(-.5*l2,.5*l2,.5*l3,.5*l2,.5*hpore-h2,size,10,1)
    surfz(-.5*R,.5*R,.5*l0,.5*R,-.5*hpore+hmem,size,5,5)
    surfz(-.5*R,.5*R,.5*l0,.5*R,-.5*hpore,size,5,5)
    #top-right
    surfz(.5*l1,.5*l0,0.,.5*l1,.5*hpore,size,5,5)
    surfz(.5*l2,.5*l1,0.,.5*l2,.5*hpore-h1,size,5,5)
    surfz(.5*l3,.5*l2,0.,.5*l3,.5*hpore-h2,size,5,5)
    surfz(.5*l0,.5*R,0.,.5*l0,-.5*hpore+hmem,size,5,5)
    surfz(.5*l0,.5*R,0.,.5*l0,-.5*hpore,size,5,5)
    #top-left
    surfz(-.5*l1,-.5*l0,0.,.5*l1,.5*hpore,size,5,5)
    surfz(-.5*l2,-.5*l1,0.,.5*l2,.5*hpore-h1,size,5,5)
    surfz(-.5*l3,-.5*l2,0.,.5*l3,.5*hpore-h2,size,5,5)
    surfz(-.5*l0,-.5*R,0.,.5*l0,-.5*hpore+hmem,size,5,5)
    surfz(-.5*l0,-.5*R,0.,.5*l0,-.5*hpore,size,5,5)
    #right
    surfx(0.,.5*l1,.5*hpore-h1,.5*hpore,.5*l1,size,5,5)
    surfx(0.,.5*l2,.5*hpore-h2,.5*hpore-h1,.5*l2,size,5,5)
    surfx(0.,.5*l3,-.5*hpore,.5*hpore-h2,.5*l3,size,5,5)
    surfx(0.,.5*l0,-.5*hpore+hmem,.5*hpore,.5*l0,size,5,5)
    #left
    surfx(0.,.5*l1,.5*hpore-h1,.5*hpore,-.5*l1,size,5,5)
    surfx(0.,.5*l2,.5*hpore-h2,.5*hpore-h1,-.5*l2,size,5,5)
    surfx(0.,.5*l3,-.5*hpore,.5*hpore-h2,-.5*l3,size,5,5)
    surfx(0.,.5*l0,-.5*hpore+hmem,.5*hpore,-.5*l0,size,5,5)


    ax.scatter(X,Y,Z)
    ax.scatter(X,-Y,Z)
    ax.scatter(-X,Y,Z)
    ax.scatter(-X,-Y,Z)
    ax.scatter(Y,X,Z)
    ax.scatter(Y,-X,Z)
    ax.scatter(-Y,X,Z)
    ax.scatter(-Y,-X,Z)
    plt.tight_layout()
    plt.show()


    plt.scatter(X_points,Y_points)
    plt.plot([0.,1.,1.,0.,0.,1.],[0.,0.,1.,1.,0.,1.],color='blue')
    ax=plt.gca()
    ax.set_aspect(1)
    plt.tight_layout()
    plt.show()
