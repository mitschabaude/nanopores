from dolfin import *
from nanopores.tools import fields
import sys
from random import gauss
import math
import numpy as np
import nanopores as nano
import nanopores.geometries.pughpore as pughpore
fields.set_dir("/tmp/nanopores/")

functions, mesh = fields.get_functions("pugh_distance")
u1 = functions["pugh_distance"]
plot(u1, interactive=True)




geop = nano.Params(pughpore.params)
physp = nano.Physics(name="pore_mol")


kT = physp.kT
eta = physp.eta
rMolecule = geop.rMolecule

Dmol = kT/(6.*math.pi*eta*rMolecule*1e-9) # [m^2/s]
gamma = (6.*math.pi*eta*rMolecule) #friction [microgramm/s]
steps = 1e8 # [ns]
tau = .05 # [ns]
C = tau/gamma*1e9 # [s^2/kg * 1e9 nm/m]
coeff = math.sqrt(2*Dmol*1e9*tau) # [nm]

F=[0.,0.,-1e-11]


X = np.zeros(steps)
Y = np.zeros(steps)
Z = np.zeros(steps)
Z[0] = 30.
