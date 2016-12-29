#from matplotlib import cm
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from nanopores.tools import fields
#import nanopores.geometries.pughpore as pughpore
import folders
from mirror import xf, Felx, Fely, Felz, Fdragx, Fdragy, Fdragz, Fx, Fy, Fz, Jf, Ja
#import nanopores

#up = nanopores.user_params(pughpore.params, k=3)


#R = up.R
#H = up.H
#l0 = up.l0
#l1 = up.l1
#l2 = up.l2
#l3 = up.l3
#l4 = up.l4
#hpore = up.hpore
#hmem = up.hmem
#h2 = up.h2
#h1 = up.h1
#h4 = up.h4
#rMolecule = up.rMolecule
#eps = 0.1
#r = rMolecule + eps
#
#fac = np.array([.5*l0*1.2,.5*l0,.5*l1-r,.5*l1-r,
#                .5*l2-r,.5*l3-r,.5*l3-r,.5*l3-r,.5*l3-r,.5*l3-r])
#heights = np.array([.5*hpore+5.,.5*hpore+rMolecule,.5*hpore,.5*(hpore-h1),
#              .5*hpore-h1,.5*hpore-h2,-.5*hpore+.75*(hpore-h2),
#              -.5*hpore+.5*(hpore-h2),-.5*hpore+.25*(hpore-h2),-.5*hpore])
#height=heights[6]

#xzi=[i for i,v in enumerate(xf) if v[2]==height]
#xz=[xf[i] for i in xzi]
#testx=[xz[i][0] for i in range(len(xz))]
#testy=[xz[i][1] for i in range(len(xz))]
#Jheight=[Jf[i] for i in xzi]
#
#fig=plt.figure()
#ax=fig.add_subplot(111,projection='3d')
#ax.scatter(testx,testy,Jheight)
#
#plt.show()
lenx=len(xf)
p=4.
maxdist=5.
def d(x,y):
    return (x[0]-y[0])**2+(x[1]-y[1])**2+(x[2]-y[2])**2
def w(i,x):
    dist=d(x,xf[i])**.5
    if dist<=maxdist:
        return 1./(d(x,xf[i])**p)
    else: return 0.
def Point(x):
    if x in xf:
        return xf.index(x)
    else: return False
def Fel(x):
    point=Point(x)
    if point is not False:
        return Felf[point]
    else:
        W = np.array([w(i,x) for i in range(lenx)])
        return [np.inner(W,Felx)/np.sum(W),np.inner(W,Fely)/np.sum(W),np.inner(W,Felz)/np.sum(W)]
def Fdrag(x):
    point=Point(x)
    if point is not False:
        return Fdragf[point]
    else:
        W = np.array([w(i,x) for i in range(lenx)])
        return [np.inner(W,Fdragx)/np.sum(W),np.inner(W,Fdragy)/np.sum(W),np.inner(W,Fdragz)/np.sum(W)]
def F(x):
    point=Point(x)
    if point is not False:
        return Ff[point]
    else:
        W = np.array([w(i,x) for i in range(lenx)])
        return [np.inner(W,Fx)/np.sum(W),np.inner(W,Fy)/np.sum(W),np.inner(W,Fz)/np.sum(W)]
def J(x):
    point=Point(x)
    if point is not False:
        return Jf[point]
    else:
        W = np.array([w(i,x) for i in range(lenx)])
        return np.inner(W,Ja)/np.sum(W)
def Fuz(x):
    point=Point(x)
    if point is not False:
        return Fz[point]
    else:
        W = np.array([w(i,x) for i in range(lenx)])
        return np.inner(W,Fz)/np.sum(W)
