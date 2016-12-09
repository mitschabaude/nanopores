#from matplotlib import cm
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from nanopores.tools import fields
#import nanopores.geometries.pughpore as pughpore
import folders
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


params=dict(bulkcon=1000.)
data=fields.get_fields("pugh",**params)
x=data["x"]
Fel=data["Fel"]
Fdrag=data["Fdrag"]
J=data["J"]
F=data["F"]

xq=list(x)
Felq=list(Fel)
Fdragq=list(Fdrag)
Jq=list(J)
Fq=list(F)
# diagonal mirror
for i in range(len(x)):
    if x[i][0]!=x[i][1]:
        xq.append([x[i][1],x[i][0],x[i][2]])
        Felq.append([Fel[i][1],Fel[i][0],Fel[i][2]])
        Fdragq.append([Fdrag[i][1],Fdrag[i][0],Fdrag[i][2]])
        Fq.append([F[i][1],F[i][0],F[i][2]])
        Jq.append(J[i])
xh=list(xq)
Felh=list(Felq)
Fdragh=list(Fdragq)
Jh=list(Jq)
Fh=list(Fq)
# left-right mirror
for i in range(len(xq)):
    if not (xq[i][0]==0. and xq[i][1]==0.):
        xh.append([-xq[i][0],xq[i][1],xq[i][2]])
        Felh.append([-Felq[i][0],Felq[i][1],Felq[i][2]])
        Fdragh.append([-Fdragq[i][0],Fdragq[i][1],Fdragq[i][2]])
        Fh.append([-Fq[i][0],Fq[i][1],Fq[i][2]])
        Jh.append(Jq[i])
xf=list(xh)
Felf=list(Felh)
Fdragf=list(Fdragh)
Jf=list(Jh)
Ff=list(Fh)
# top-bottom mirror
for i in range(len(xh)):
    if not (xh[i][0]==0. and xh[i][1]==0.):
        xf.append([xh[i][0],-xh[i][1],xh[i][2]])
        Felf.append([Felh[i][0],-Felh[i][1],Felh[i][2]])
        Fdragf.append([Fdragh[i][0],-Fdragh[i][1],Fdragh[i][2]])
        Ff.append([Fh[i][0],-Fh[i][1],Fh[i][2]])
        Jf.append(Jh[i])
lenx=len(xf)

Felx=np.array([Felf[i][0] for i in range(lenx)])
Fely=np.array([Felf[i][1] for i in range(lenx)])
Felz=np.array([Felf[i][2] for i in range(lenx)])
Fdragx=np.array([Fdragf[i][0] for i in range(lenx)])
Fdragy=np.array([Fdragf[i][1] for i in range(lenx)])
Fdragz=np.array([Fdragf[i][2] for i in range(lenx)])
Fx=np.array([Ff[i][0] for i in range(lenx)])
Fy=np.array([Ff[i][1] for i in range(lenx)])
Fz=np.array([Ff[i][2] for i in range(lenx)])
Ja=np.array(Jf)
