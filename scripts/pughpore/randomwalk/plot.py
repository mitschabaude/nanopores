import nanopores as nano
import nanopores.geometries.pughpore as pughpore
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
#from numpy.random import random
import matplotlib.pyplot as plt
import nanopores.geometries.pughpore as pughpore
import nanopores
import sys
import os
import nanopores.tools.fields as f
HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")
DATADIR = os.path.join(HOME, "Dropbox", "Paper Howorka", "data", "fields")
f.set_dir(DATADIR)
import sys
if len(sys.argv)==1:
    sys.exit('Integer missing! Type:\npython plot.py i\nwith i=integer.')
i=int(sys.argv[1])

up = nano.Params(pughpore.params, k=3)
hpore=up.hpore

#X = np.load('X.npy')
#Y = np.load('Y.npy')
#Z = np.load('Z.npy')
#J = np.load('J.npy')
#T = np.load('T.npy')
params=dict(avgbind=1e7,P_bind=3.e-4,z0=hpore/2.+5.)
data=f.get_fields("randomwalk2",**params)
X = np.array(data["X"][i])
Y = np.array(data["Y"][i])
Z = np.array(data["Z"][i])
T = np.array(data["T"][i])
J = np.array(data["J"][i])
long = np.where(T>100.)
amplitude = 2060.-np.inner(J,T)/np.sum(T)
for i in range(1,T.shape[0]):
    T[i]=T[i]+T[i-1]
tau_off=T[-1]
print 'tau_off = %.3f ms'% (tau_off*1e-6)
print 'amplitude = %.0f pA'% amplitude
#J_a = J[0]
#J_b = J[-1]
J=np.append(np.array([2060.,2060.]),J)
J=np.append(J,np.array([2060.,2060.]))
T=np.append(np.array([-1e9,0.]),T)
T=np.append(T,np.array([tau_off,1e9+tau_off]))
T=T*1e-9

#R = up.R
R = 50.
#H = up.H
H = 100.
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

def surfx(y1,y2,z1,z2,d,size,rs,cs):
    Y = np.linspace(y1,y2,size)
    Z = np.linspace(z1,z2,size)
    Y, Z = np.meshgrid(Y,Z)
    X = np.zeros(size)+d
    surf = ax.plot_surface(X,Y,Z, rstride=rs, cstride=cs, alpha=alpha,color=color)

def surfy(x1,x2,z1,z2,d,size,rs,cs):
    X = np.linspace(x1,x2,size)
    Z = np.linspace(z1,z2,size)
    X, Z = np.meshgrid(X,Z)
    Y = np.zeros(size)+d
    surf = ax.plot_surface(X,Y,Z, rstride=rs, cstride=cs, alpha=alpha,color=color)

def surfz(x1,x2,y1,y2,d,size,rs,cs):
    X = np.linspace(x1,x2,size)
    Y = np.linspace(y1,y2,size)
    X, Y = np.meshgrid(X,Y)
    Z = np.zeros(size)+d
    surf = ax.plot_surface(X,Y,Z, rstride=rs, cstride=cs, alpha=alpha,color=color)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.view_init(elev=0,azim=270)
#ax.set_xlim([-.5*R,.5*R])
#ax.set_ylim([-.5*R,.5*R])
#ax.set_zlim([-.5*H,.5*H])
ax.set_aspect(3.)

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

plt.plot(X,Y,Z)
ax.scatter(X[long],Y[long],Z[long],c='r',s=50)

plt.tight_layout()
plt.show()
plt.plot(T,J)
ax=plt.gca()
ax.set_xlabel('Time [s]')
ax.set_ylabel('Current [pA]')
ax.set_ylim([1950,2100])
plt.tight_layout()
plt.show()
