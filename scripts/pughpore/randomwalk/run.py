from dolfin import *
from nanopores.tools import fields
import sys
from random import gauss
import math
import numpy as np
import nanopores as nano
import nanopores.geometries.pughpore as pughpore
fields.set_dir("/home/benjamin/Desktop/data/")
from smooth_D import f, fp

functions, mesh = fields.get_functions("pugh_distance")
dis = functions["pugh_distance"]
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
steps = 1e5 # [ns]
tau = .05 # [ns]
C = tau/gamma*1e9 # [s^2/kg * 1e9 nm/m]
coeff = math.sqrt(2*Dmol*1e9*tau) # [nm]

F=[0.,0.,-1e-11]
#F=[0.,0.,0.]


X = np.zeros(steps)
Y = np.zeros(steps)
Z = np.zeros(steps)
Z[0] = 10.
i=0
bond = 0
while i<steps-1:
    xi_x=gauss(0.,1.)
    xi_y=gauss(0.,1.)
    xi_z=gauss(0.,1.)
    Force = F
    Dfac = f(Z[i])
    X[i+1] = X[i] + coeff*xi_x + C*Force[0]
    Y[i+1] = Y[i] + coeff*xi_y + C*Force[1]
    Z[i+1] = Z[i] + coeff*xi_z*math.sqrt(Dfac) + C*Force[2]*Dfac + fp(Z[i])
    if dis(argument(X[i+1],Y[i+1],Z[i+1])) < rMolecule:
        i-=1
        bond+=1
    i+=1
    if not (Z[i]<=H/2. and Z[i]>=-H/2 and X[i] <=R/2 and X[i] >=-R/2 and Y[i] <=R/2 and Y[i] >=-R/2):
        X = X[:i+1]
        Y = Y[:i+1]
        Z = Z[:i+1]
        break
print 'bindings',bond
np.save('X',X)
np.save('Y',Y)
np.save('Z',Z)
