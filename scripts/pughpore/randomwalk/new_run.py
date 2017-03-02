import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from dolfin import *
from nanopores.tools import fields
import sys
from random import gauss, expovariate
import math
from math import atan, pi, atan2, sqrt
import numpy as np
import nanopores as nano
import nanopores.geometries.pughpore as pughpore
from get_F import Force, Current
from get_D import Dx, Dy, Dz, dxDx, dyDy, dzDz, dis
import os
from time import time as timer
sys.exit()


HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")

DATADIR = os.path.join(HOME,"Dropbox", "nanopores", "fields")

import nanopores.tools.fields as fields
fields.set_dir(DATADIR)

def argument(x,y,z):
    return np.array([float(x),float(y),float(z)])

geop = nano.Params(pughpore.params)
physp = nano.Physics(name="pore_mol")

kT = physp.kT
eta = physp.eta
#H = geop.H
#R = geop.R
H = 100.
R = 50.


l0 =        geop.l0
l1 =        geop.l1
l2 =        geop.l2
l3 =        geop.l3
l4 =        geop.l4
hpore =     geop.hpore
hmem =      geop.hmem
h2 =        geop.h2
h1 =        geop.h1
h4 =        geop.h4
rMolecule = geop.rMolecule
eps = 0.1
beps = (l3 - rMolecule)*1e-1
r = rMolecule + eps
p0=hpore/2.
p1=p0-h1
p2=p0-h2
p3=-hpore/2.


def R_(z):
    if z>=p3 and z<=p2:
        return l3/2.
    elif z>=p2 and z<=p1:
        return l2/2.
    elif z>=p1 and z<=p0:
        return l0/2.
    else: return R/2.

#def fac(z):
#    if z>=p3 and z<=p2:
#        return l3/2.-r
#    elif z>p2 and z<p2+r:
#        x=z-p2
#        return -sqrt(r**2-x**2)+l3/2.
#    elif z>=p2+r and z<=p1:
#        return l2/2.-r
#    elif z>p1 and z<p1+r:
#        x=z-p1
#        return -sqrt(r**2-x**2)+l2/2.
#    elif z>=p1+r and z<=p0:
#        return l0/2.-r
#    elif z>p0 and z<p0+r:
#        x=z-p0
#        return -sqrt(r**2-x**2)+l0/2.
#    elif z<p3 and z>p3-r:
#        x=z-p3
#        return -sqrt(r**2-x**2)+l3/2.
#    else: return R/2.

#params=dict(avgbind1=7e6,avgbind2=1e2,P_bind1=5.e-3,P_bind2=1e-1,z0=hpore/2.+5.)

Dmol = kT/(6.*math.pi*eta*rMolecule*1e-9) # [m^2/s]
gamma = (6.*math.pi*eta*rMolecule) #friction [microgramm/s]
maxiter = 1e6 # [ns]
tau = .05 # [ns]
C = tau/gamma*1e9 # [s^2/kg * 1e9 nm/m]
coeff = math.sqrt(2*Dmol*1e9*tau) # [nm]
#avgbinding = 10000000.
#P_bind = 3.e-4

F=[0.,0.,-1e-11]
#F=[0.,0.,0.]

def hatfct(ang):
    x=(ang+2*pi)%(pi/2.)
    if x<=pi/4.:
        return x
    else:
        return pi/2.-x
def D(x,y,z):
    if z>hpore/2. or z<-hpore/2.:
        return [[1.,1.,1.],[0.,0.,0.]]
    else:
        if x==0 and y==0:
            return [[Dx(0.),Dy(0.),Dz(0.)],[dDx(0.),dDy(0.),dDz(0.)]]
        else:
            ang=atan2(y,x)
            ang2=hatfct(ang)
            A=np.array([[cos(ang),-sin(ang),0.],[sin(ang),cos(ang),0.],[0.,0.,1.]])
            dist=sqrt(x**2+y**2)*cos(ang2)/(R_(z))
            vec1=A.dot(np.array([Dx(dist),Dy(dist),Dz(dist)]))
            vec2=A.dot(np.array([dDx(dist),dDy(dist),dDz(dist)]))
            return [list(vec1),list(vec2)]

def run(params):
    z0 = params["z0"]
    X = np.array([0.])
    Y = np.array([0.])
    Z = np.array([z0])
    J1 = np.array([])
    T = np.array([])
    bind1 = 0
    bind2 = 0
    avgbind1=params["avgbind1"]
    P_bind1=params["P_bind1"]
    avgbind2=params["avgbind2"]
    P_bind2=params["P_bind2"]
    ffa = True
    i=0
    ood = False
    while i<maxiter and Z[-1]>=-hpore/2.-2.:
        if ood:
	    bind1 = 0
	    bind2 = 0
            i=0
            ood = False
            ffa = True
            X = np.array([0.])
            Y = np.array([0.])
            Z = np.array([z0])
            T = np.array([])
            J1 = np.array([])
        add=tau
        xi_x=gauss(0.,1.)
        xi_y=gauss(0.,1.)
        xi_z=gauss(0.,1.)
        Force = F
#	[[Dxfac, Dyfac, Dzfac],[DDx,DDy,DDz]]=D(X[-1],Y[-1],Z[-1])
#        x_new = X[-1] + coeff*xi_x*math.sqrt(abs(Dxfac)) + C*Force[0]*Dxfac + DDx*tau*Dmol
#        y_new = Y[-1] + coeff*xi_y*math.sqrt(abs(Dyfac)) + C*Force[1]*Dyfac + DDy*tau*Dmol
#        z_new = Z[-1] + coeff*xi_z*math.sqrt(abs(Dzfac)) + C*Force[2]*Dzfac + DDz*tau*Dmol
        x_new = X[-1] + coeff*xi_x + C*Force[0]
        y_new = Y[-1] + coeff*xi_y + C*Force[1]
        z_new = Z[-1] + coeff*xi_z + C*Force[2]
        if dis(argument(x_new,y_new,z_new)) < rMolecule:
            x_new = X[-1]
            y_new = Y[-1]
            z_new = Z[-1]
            if ffa and np.random.binomial(1,P_bind1)==1 and Z[-1]<=hpore/2.-h2-5 and Z[-1]>=-hpore/2.+h4:
                add+=expovariate(lambd=1./avgbind1)
#                print add
                bind1+=1
            elif ffa and np.random.binomial(1,P_bind2)==1 and ((Z[-1]<=-hpore/2.+h4 and Z[-1]>=-hpore/2.+0.) or (Z[-1]<=hpore/2.-h2 and Z[-1]>=hpore/2.-h2-5.)):
                add+=expovariate(lambd=1./avgbind2)
                bind2+=1
            else:
                add+=0.
            ffa = False
        elif dis(argument(x_new,y_new,z_new)) < rMolecule + beps:
            pass
        else:
            ffa = True
        X = np.append(X,x_new)
        Y = np.append(Y,y_new)
        Z = np.append(Z,z_new)
        if abs(Z[-1])>30.:
            print 'Out of domain!'
            ood = True
            Z[-1]=0.
        J1=np.append(J1,J(Z[-1]))
        T =np.append(T,add)
        i+=1
        if not (Z[i]<=H/2. and Z[i]>=-H/2 and X[i] <=R/2 and X[i] >=-R/2 and Y[i] <=R/2 and Y[i] >=-R/2):
            break
    if i>=maxiter:
        print 'randomwalk: more than 1e6 steps!'
#    X=[list(X)]
#    Y=[list(Y)]
#    Z=[list(Z)]
#    T=[list(T)]
#    J1=[list(J1)]
    tau_off = np.sum(T)*1e-6
    amp = (2060.-np.inner(T*1e-6,J1)/tau_off)/2060.*100.
    if tau_off<1.:
        t1 = [tau_off]
        a1 = [amp]
        t2 = []
        a2 = []
    else:
        t2 = [tau_off]
        a2 = [amp]
        t1 = []
        a1 = []
        
#    print 'bind1 = %i'%bind1
#    print 'bind2 = %i'%bind2
#    print '_'
    fields.save_fields("rw_2_bind",params,t1=t1,a1=a1,t2=t2,a2=a2)
