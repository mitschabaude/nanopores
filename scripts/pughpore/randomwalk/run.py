import matplotlib.pyplot as plt
from dolfin import *
from nanopores.tools import fields
import sys
from random import gauss, expovariate
import math
import numpy as np
import nanopores as nano
import nanopores.geometries.pughpore as pughpore
from get_D import f, fp
from get_J import Jf as J
import os
from time import time as timer


HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")

DATADIR = os.path.join(HOME, "Dropbox", "Paper Howorka", "data", "fields")

import nanopores.tools.fields as fields
fields.set_dir(DATADIR)

functions, mesh = fields.get_functions("pugh_distance", h=.5)
dis = functions["y"]
def argument(x,y,z):
    return np.array([float(x),float(y),float(z)])

geop = nano.Params(pughpore.params)
physp = nano.Physics(name="pore_mol")

kT = physp.kT
eta = physp.eta
rMolecule = geop.rMolecule
#H = geop.H
#R = geop.R
H = 100.
R = 50.
hpore = geop.hpore
h2 = geop.h2

params=dict(avgbind=1e7,P_bind=3.e-4,z0=hpore/2.+5.)

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


def run(params=params):
    X = np.array([0.])
    Y = np.array([0.])
    Z = np.array([params["z0"]])
    J1 = np.array([])
    T = np.array([])
    avgbind=params["avgbind"]
    P_bind=params["P_bind"]
    i=0
    while i<maxiter and Z[-1]>=-hpore/2.-2.:
        xi_x=gauss(0.,1.)
        xi_y=gauss(0.,1.)
        xi_z=gauss(0.,1.)
        Force = F
        Dfac = f(Z[-1])
        x_new = X[-1] + coeff*xi_x + C*Force[0]
        y_new = Y[-1] + coeff*xi_y + C*Force[1]
        z_new = Z[-1] + coeff*xi_z*math.sqrt(Dfac) + C*Force[2]*Dfac + fp(Z[i])*Dmol*tau
        if dis(argument(x_new,y_new,z_new)) < rMolecule:
            x_new = X[-1]
            y_new = Y[-1]
            z_new = Z[-1]
            if np.random.binomial(1,P_bind)==1 and Z[-1]<=hpore/2.-h2 and Z[-1]>=-hpore/2.+1.:
                add=expovariate(lambd=1./avgbind)
            else:
                add=0.
        else:
            add=tau
        X = np.append(X,x_new)
        Y = np.append(Y,y_new)
        Z = np.append(Z,z_new)
        J1=np.append(J1,J(Z[-1]))
        T =np.append(T,add)
        i+=1
        if not (Z[i]<=H/2. and Z[i]>=-H/2 and X[i] <=R/2 and X[i] >=-R/2 and Y[i] <=R/2 and Y[i] >=-R/2):
            break
    if i>=maxiter:
        print 'randomwalk: more than 1e6 steps!'
    X=[list(X)]
    Y=[list(Y)]
    Z=[list(Z)]
    T=[list(T)]
    J1=[list(J1)]
    fields.save_fields("randomwalk3",params,X=X,Y=Y,Z=Z,T=T,J=J1)
