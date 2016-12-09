from matplotlib.collections import PolyCollection
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import nanopores.geometries.pughpore as pughpore
import nanopores
from inv_interp import Fel, Fdrag, F, J, Fuz

fig = plt.figure()
ax = fig.gca(projection='3d')

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

verts=[]
xs = [l3*.5,R*.5,R*.5,l4*.5,l4*.5,l0*.5,l0*.5,l1*.5,l1*.5,l2*.5,l2*.5,l3*.5]
xsm = [-xs[i] for i in range(len(xs))]
ys = [-hpore*.5,-hpore*.5,-hpore*.5+hmem,-hpore*.5+hmem,-hpore*.5+h4,-hpore*.5+h4,hpore*.5,hpore*.5,hpore*.5-h1,hpore*.5-h1,hpore*.5-h2,hpore*.5-h2]
verts.append(list(zip(xs,ys)))
verts.append(list(zip(xsm,ys)))


poly = PolyCollection(verts)
poly.set_alpha(.3)

ax.set_xlim3d(-15, 15)
ax.set_ylim3d(-15, 15)
ax.set_zlim3d(0.,2e-9)
lx=40*2
ly=80*2
X_=np.linspace(-R*.5,R*.5,lx)
Y_=np.linspace(-hpore*.5,hpore*.5,ly)
X, Y = np.meshgrid(X_,Y_)
#np.save('X',X)
#np.save('Y',Y)
#ZJ=np.array([[Fuz([X_[i],0.,Y_[j]]) for i in range(X_.shape[0])] for j in range(Y_.shape[0])])
#np.save('Fz2',ZJ)
ZJ=np.load('ZJ4.npy')
surf=ax.plot_surface(X,Y,ZJ, rstride=1, cstride=1, cmap=cm.viridis, linewidth=0., alpha=1.)
fig.colorbar(surf, shrink=0.5, aspect=5)

ax.add_collection3d(poly, zs=[0.,0.], zdir='z')

#ax.view_init(23.,-32.)
ax.view_init(90.,-90.)
plt.tight_layout()
plt.show()
