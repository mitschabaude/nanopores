from math import isnan
import numpy as np
try:
    xf = np.load('xf.npy')
    Fx = np.load('Fx.npy')
    Fy = np.load('Fy.npy')
    Fz = np.load('Fz.npy')
    Ja = np.load('Ja.npy')
except:
    import os
    import nanopores.tools.fields as fields
    HOME = os.path.expanduser("~")
    PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
    FIGDIR = os.path.join(PAPERDIR, "figures", "")

    DATADIR = os.path.join(HOME,"Dropbox", "nanopores", "fields")

    fields.set_dir(DATADIR)

    params=dict(Nmax=500000.)
    data = fields.get_fields("pugh_force",**params)

    x=data["x"]
    #Fel=data["Fel"]
    #Fdrag=data["Fdrag"]
    J=data["J"]
    F=data["F"]

    xq=list(x)
    #Felq=list(Fel)
    #Fdragq=list(Fdrag)
    Jq=list(J)
    Fq=list(F)
    # diagonal mirror
    for i in range(len(x)):
        if x[i][0]!=x[i][1]:
            xq.append([x[i][1],x[i][0],x[i][2]])
    #        Felq.append([Fel[i][1],Fel[i][0],Fel[i][2]])
    #        Fdragq.append([Fdrag[i][1],Fdrag[i][0],Fdrag[i][2]])
            Fq.append([F[i][1],F[i][0],F[i][2]])
            Jq.append(J[i])
    xh=list(xq)
    #Felh=list(Felq)
    #Fdragh=list(Fdragq)
    Jh=list(Jq)
    Fh=list(Fq)
    # left-right mirror
    for i in range(len(xq)):
        if not (xq[i][0]==0. and xq[i][1]==0.):
            xh.append([-xq[i][0],xq[i][1],xq[i][2]])
    #        Felh.append([-Felq[i][0],Felq[i][1],Felq[i][2]])
    #        Fdragh.append([-Fdragq[i][0],Fdragq[i][1],Fdragq[i][2]])
            Fh.append([-Fq[i][0],Fq[i][1],Fq[i][2]])
            Jh.append(Jq[i])
    xf=list(xh)
    #Felf=list(Felh)
    #Fdragf=list(Fdragh)
    Jf=list(Jh)
    Ff=list(Fh)
    # top-bottom mirror
    for i in range(len(xh)):
        if not (xh[i][0]==0. and xh[i][1]==0.):
            xf.append([xh[i][0],-xh[i][1],xh[i][2]])
    #        Felf.append([Felh[i][0],-Felh[i][1],Felh[i][2]])
    #        Fdragf.append([Fdragh[i][0],-Fdragh[i][1],Fdragh[i][2]])
            Ff.append([Fh[i][0],-Fh[i][1],Fh[i][2]])
            Jf.append(Jh[i])
    lenx=len(xf)

    #Felx=np.array([Felf[i][0] for i in range(lenx)])
    #Fely=np.array([Felf[i][1] for i in range(lenx)])
    #Felz=np.array([Felf[i][2] for i in range(lenx)])
    #Fdragx=np.array([Fdragf[i][0] for i in range(lenx)])
    #Fdragy=np.array([Fdragf[i][1] for i in range(lenx)])
    #Fdragz=np.array([Fdragf[i][2] for i in range(lenx)])
    xf=np.array(xf)
    Fx=np.array([Ff[i][0] for i in range(lenx)])
    Fy=np.array([Ff[i][1] for i in range(lenx)])
    Fz=np.array([Ff[i][2] for i in range(lenx)])
    Ja=np.array(Jf)
    np.save('xf',xf)
    np.save('Fx',Fx)
    np.save('Fy',Fy)
    np.save('Fz',Fz)
    np.save('Ja',Ja)

from  scipy.interpolate import LinearNDInterpolator
Fxi = LinearNDInterpolator(xf,Fx)
Fyi = LinearNDInterpolator(xf,Fy)
Fzi = LinearNDInterpolator(xf,Fz)
def Force(x,y,z):
    ret = [Fxi.__call__(np.array([x,y,z]))[0],Fyi.__call__(np.array([x,y,z]))[0],Fzi.__call__(np.array([x,y,z]))[0]]
    if isnan(ret[0]) or isnan(ret[1]) or isnan(ret[2]):
        return [0.,0.,0.]
    else: return ret
J = LinearNDInterpolator(xf,Ja)
def Current(x,y,z):
    if z>21. or z<-23-1.6*2.0779:
        return 7.523849e-10
    elif z<-23.:
        k = (6.186351-7.523849)*1e-10/(1.6*2.0779)
        d = 6.186351*1e-10-k*(-23.)
        return k*z+d
    else:
        return -J.__call__(np.array([x,y,z]))[0]


#fig = plt.figure()
#ax = fig.gca(projection='3d')
#
#up = nanopores.user_params(pughpore.params, k=3)
#
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
#verts=[]
#xs = [l3*.5,R*.5,R*.5,l4*.5,l4*.5,l0*.5,l0*.5,l1*.5,l1*.5,l2*.5,l2*.5,l3*.5]
#xsm = [-xs[i] for i in range(len(xs))]
#ys = [-hpore*.5,-hpore*.5,-hpore*.5+hmem,-hpore*.5+hmem,-hpore*.5+h4,-hpore*.5+h4,hpore*.5,hpore*.5,hpore*.5-h1,hpore*.5-h1,hpore*.5-h2,hpore*.5-h2]
#verts.append(list(zip(xs,ys)))
#verts.append(list(zip(xsm,ys)))
#
#
#poly = PolyCollection(verts)
#poly.set_alpha(.3)
#
#ax.set_xlim3d(-15, 15)
#ax.set_ylim3d(-15, 15)
#ax.set_zlim3d(1.1e-9,1.3e-9)
#lx=80
#ly=160
#X_=np.linspace(-R*.5,R*.5,lx)
#Y_=np.linspace(-hpore*.5,hpore*.5,ly)
##X_=np.linspace(-(l3*.5-r),l3*.5-r)
#X, Y = np.meshgrid(X_,Y_)
##np.save('X',X)
##np.save('Y',Y)
##ZJ=np.array([[f.__call__(np.array([X_[i],0.,Y_[j]]))[0] for i in range(X_.shape[0])] for j in range(Y_.shape[0])])
##np.save('J_int_shrink',ZJ)
#ZJ=np.load('J_int.npy')
#surf=ax.plot_surface(X,Y,ZJ, rstride=1, cstride=1, cmap=cm.viridis, linewidth=0., alpha=1.)
#fig.colorbar(surf, shrink=0.5, aspect=5)
#
#ax.add_collection3d(poly, zs=[0.,0.], zdir='z')
#
#ax.view_init(23.,-32.)
##ax.view_init(90.,-90.)
#plt.tight_layout()
#plt.show()
