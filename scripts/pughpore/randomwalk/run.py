import sys
from random import gauss
import math
import numpy as np
import nanopores as nano
import nanopores.geometries.pughpore as pughpore

geop = nano.Params(pughpore.params)
physp = nano.Physics(name="pore_mol")

kT = physp.kT
eta = physp.eta
rMolecule = geop.rMolecule

Dmol = kT/(6.*math.pi*eta*rMolecule*1e-9)





###########################################################
#steps=
#X = np.zeros(steps)
#Y = np.zeros(steps)
#Z = np.zeros(steps)
