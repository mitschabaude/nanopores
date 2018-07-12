from matplotlib import gridspec
import matplotlib.path as mplPath
import numpy as np
from get_F import Force, Current
from nanopores.models.pughpore import polygon
from nanopores.models.pughpoints import plot_polygon
import matplotlib.pyplot as plt
import nanopores as nano
import nanopores.geometries.pughpore as pughpore
#import os

up    = nano.Params(pughpore.params, k=3)

l0    = up.l0
l1    = up.l1
l2    = up.l2
l3    = up.l3
h1    = up.h1
h2    = up.h2
hpore = up.hpore

#import nanopores.tools.fields as f
#HOME = os.path.expanduser("~")
#PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
#FIGDIR = os.path.join(PAPERDIR, "figures", "")
#DATADIR = os.path.join(HOME,"Dropbox", "nanopores", "fields")
#f.set_dir(DATADIR)


#    30
#
#
#-10     10
#
#   -25

extent=[-10,10,-25,30]
vmin = 5.9109868319135862e-10
vmax = 7.5238490173716702e-10
poly=polygon(rmem=20.)
Polygon=mplPath.Path(np.array([[poly[i][0],poly[i][1]] for i in range(len(poly))]))

try:
    matrix1=np.load('matrix1.npy')
    matrix2=np.load('matrix2.npy')
except:
    hinv=5
    xlen=(extent[1]-extent[0])*hinv
    ylen=(extent[3]-extent[2])*hinv
    matrix1=np.zeros((ylen,xlen))
    matrix2=np.zeros((ylen,xlen))
    for i in range(ylen):
        for j in range(xlen):
            x = float(j)/hinv+extent[0]
            y = extent[3]-float(i)/hinv
            if Polygon.contains_point((x,y)) or Polygon.contains_point((-x,y)):# or (y<-10. and abs(x)>=7.):
                matrix1[i][j]=float('nan')
                matrix2[i][j]=float('nan')
            else:
                matrix1[i][j]=Current(x,0.,y)
                matrix2[i][j]=(vmax-Current(x,0.,y))/vmax
    np.save('matrix1',matrix1)
    np.save('matrix2',matrix2)


#args=dict()
args=dict(extent=extent)
#args=dict(extent=extent,vmin=vmin,vmax=vmax)

fig=plt.figure(figsize=(20,10),dpi=80)
gs=gridspec.GridSpec(1, 2, width_ratios=[1,1])

ax0=plt.subplot(gs[0])


plt.title('Current')
ax=plt.gca()
ax.set_aspect('equal')
ax.set_xlim([-20.,20.])
ax.set_ylim([-23.,24.])
ax.set_xticks(np.arange(-20.,20.,2.))
ax.set_yticks(np.arange(-23.,24.,2.))
ax.grid()
data = ax0.imshow(matrix1,cmap=plt.cm.viridis,interpolation='none',**args)
cbar = fig.colorbar(data)
plot_polygon(ax,poly)


ax1=plt.subplot(gs[1])
plt.title('Drop')
ax=plt.gca()
ax.set_aspect('equal')
ax.set_xlim([-20.,20.])
ax.set_ylim([-23.,24.])
ax.set_xticks(np.arange(-20.,20.,2.))
ax.set_yticks(np.arange(-23.,24.,2.))
ax.grid()

data = ax1.imshow(matrix2,cmap=plt.cm.viridis,interpolation='none',**args)
cbar = fig.colorbar(data)



plot_polygon(ax,poly)
plt.tight_layout()
plt.show()
