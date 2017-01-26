from matplotlib import pyplot as plt
import numpy as np
import nanopores
import os
from nanopores.tools import fields
from scipy.interpolate import interp1d

HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")

DATADIR = os.path.join(HOME, "Dropbox", "nanopores", "fields")

fields.set_dir(DATADIR)

data = fields.get_fields("pugh_diff3D_cross", bulkbc=True, rMolecule=2.0779)
def smooth3(l):
    A=np.array(l)
    B=A[:]
    ker=np.array([1./3,1./3,1./3])
    n=int(ker.shape[0]/2.)
    for i in range(n,A.shape[0]-n):
        B[i]=np.inner(A[i-n:i+n+1],ker)
    return list(B)
def smooth5(l):
    A=np.array(l)
    B=A[:]
    ker=np.array([.2,.2,.2,.2,.2])
    n=int(ker.shape[0]/2.)
    for i in range(n,A.shape[0]-n):
        B[i]=np.inner(A[i-n:i+n+1],ker)
    return list(B)
def smootha(l):
    A=np.array(l)
    B=A[:]
    ker=np.array([10.,12.,15.,12.,10.])
    ker=ker/np.sum(ker)
    n=int(ker.shape[0]/2.)
    for i in range(n,A.shape[0]-n):
        B[i]=np.inner(A[i-n:i+n+1],ker)
    return list(B)

x = [z[0] for z in data["x"]]
data, x = fields._sorted(data, x)
eps=5e-3
x_=x[:]
#x_.extend([1.,1.+eps,1.+2*eps,1.+3*eps])
x.extend([(x[-1]+1.)/2.,1.,1.+eps,1.+2*eps,1.+3*eps,1.+4*eps,1.+5*eps])
dstr = ["x", "y", "z"]
Dxx = [D[0][0] for D in data["D"]]
Dyy = [D[1][1] for D in data["D"]]
Dzz = [D[2][2] for D in data["D"]]
Dxx_ = [D[0][0] for D in data["D"]]
Dyy_ = [D[1][1] for D in data["D"]]
Dzz_ = [D[2][2] for D in data["D"]]
Dxx.extend([0.,0.,0.,0.,0.,0.,0.])
Dyy.extend([Dyy[-1]/2.,0.,0.,0.,0.,0.,0.])
Dzz.extend([Dzz[-1]/2.,0.,0.,0.,0.,0.,0.])
#Dxx_.extend([0.,0.,0.,0.])
#Dyy_.extend([0.,0.,0.,0.])
#Dzz_.extend([0.,0.,0.,0.])
Dxx=smooth5(smooth3(Dxx))
Dyy=smooth5(smooth3(Dyy))
Dzz=smooth5(smooth3(Dzz))
Dx = interp1d(x,Dxx)
Dy = interp1d(x,Dyy)
Dz = interp1d(x,Dzz)

xc=np.linspace(0.,1.,100)


plt.plot(x_,Dxx_,color='blue',linestyle=':')
plt.scatter(x_,Dxx_,color='blue')
plt.scatter(x,Dxx,color='blue')
plt.plot(x,Dxx,color='blue',label=r"$D_{%s%s}$" % (dstr[0], dstr[0]))
plt.plot(xc,Dx(xc),color='blue')

DDxx = [0.]+[(Dxx[i+1]-Dxx[i-1])/(x[i+1]-x[i-1]) for i in range(1,len(x)-1)]+[0.]
plt.scatter(x,DDxx,color='blue')
plt.plot(x,DDxx,color='blue')



plt.plot(x_,Dyy_,color='red',linestyle=':')
plt.scatter(x_,Dyy_,color='red')
plt.scatter(x,Dyy,color='red')
plt.plot(x,Dyy,color='red',label=r"$D_{%s%s}$" % (dstr[1], dstr[1]))
plt.plot(xc,Dy(xc),color='red')

DDyy = [0.]+[(Dyy[i+1]-Dyy[i-1])/(x[i+1]-x[i-1]) for i in range(1,len(x)-1)]+[0.]
plt.scatter(x,DDyy,color='red')
plt.plot(x,DDyy,color='red')



plt.plot(x_,Dzz_,color='green',linestyle=':')
plt.scatter(x_,Dzz_,color='green')
plt.scatter(x,Dzz,color='green')
plt.plot(x,Dzz,color='green',label=r"$D_{%s%s}$" % (dstr[2], dstr[2]))
plt.plot(xc,Dz(xc),color='green')

DDzz = [0.]+[(Dzz[i+1]-Dzz[i-1])/(x[i+1]-x[i-1]) for i in range(1,len(x)-1)]+[0.]
plt.scatter(x,DDzz,color='green')
plt.plot(x,DDzz,color='green')



plt.xlabel('distance from pore center [nm]')
plt.ylabel('diffusivity relative to bulk')
plt.legend(loc='lower left')
plt.tight_layout()
plt.show()
