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
H = 80.
R = 30.

Dmol = kT/(6.*math.pi*eta*rMolecule*1e-9) # [m^2/s]
gamma = (6.*math.pi*eta*rMolecule) #friction [microgramm/s]
steps = 1e6 # [ns]
tau = .05 # [ns]
C = tau/gamma*1e9 # [s^2/kg * 1e9 nm/m]
coeff = math.sqrt(2*Dmol*1e9*tau) # [nm]
avgbinding = 25000000.
P_bind = 7.e-4

F=[0.,0.,-1e-11]
#F=[0.,0.,0.]


def run(plotres=False):
    print 'RUNNING'
    X = np.zeros(steps)
    Y = np.zeros(steps)
    Z = np.zeros(steps)
    J1 = np.array([])
    T = np.array([])
    Z[0] = 25.
    i=0
    bind = 0
    contact = 0
    time=0.
    while i<steps-1 and Z[i]>=-24.:
        xi_x=gauss(0.,1.)
        xi_y=gauss(0.,1.)
        xi_z=gauss(0.,1.)
        Force = F
        Dfac = f(Z[i])
        X[i+1] = X[i] + coeff*xi_x + C*Force[0]
        Y[i+1] = Y[i] + coeff*xi_y + C*Force[1]
        Z[i+1] = Z[i] + coeff*xi_z*math.sqrt(Dfac) + C*Force[2]*Dfac + fp(Z[i])*Dmol
        time+=.05
        if dis(argument(X[i+1],Y[i+1],Z[i+1])) < rMolecule:
            X[i+1] = X[i]
            Y[i+1] = Y[i]
            Z[i+1] = Z[i]
            contact+=1
            if np.random.binomial(1,P_bind)==1:
                bind+=1
                add=expovariate(lambd=1./avgbinding)
                time+=add
                J1=np.append(J1,np.array([J(Z[i])])) 
                T =np.append(T,np.array([add]))
        else:
            J1=np.append(J1,np.array([J(Z[i])]))
            T =np.append(T,np.array([tau]))
        i+=1
        if not (Z[i]<=H/2. and Z[i]>=-H/2 and X[i] <=R/2 and X[i] >=-R/2 and Y[i] <=R/2 and Y[i] >=-R/2):
            break
    X = X[:i+1]
    Y = Y[:i+1]
    Z = Z[:i+1]
    print 'contacts',contact
    print 'bindings',bind
    np.save('X',X)
    np.save('Y',Y)
    np.save('Z',Z)
    np.save('J',J1)
    np.save('T',T)
    tau_off=np.sum(T)
    amplitude = (2060.-np.inner(J1,T)/tau_off)/2060.*100
    print 'tau_off = %.3f ms'% (tau_off*1e-6)
    print 'A/I_0 = %.1f %%'% amplitude
    np.save('TAU',np.append(np.load('TAU.npy'),tau_off))
    np.save('A',np.append(np.load('A.npy'),amplitude))
    if plotres:
        for i in range(1,T.shape[0]):
            T[i]=T[i]+T[i-1]
        J1=np.append(np.array([2060.,2060.]),J1)
        J1=np.append(J1,np.array([2060.,2060.]))
        T=np.append(np.array([-1e9,0.]),T)
        T=np.append(T,np.array([tau_off,1e9+tau_off]))
        T=T*1e-9
        plt.plot(T,J1)
        ax=plt.gca()
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Current [pA]')
        ax.set_ylim([1950,2100])
        plt.tight_layout()
        plt.show()

