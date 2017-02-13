import numpy as np
import os
from math import pi, sinh, acosh
import nanopores as nano
import nanopores.geometries.pughpore as pughpore
import nanopores.tools.fields as fields
geop = nano.Params(pughpore.params)
physp = nano.Physics(name="pore_mol")

HOME = os.path.expanduser("~")
PAPERDIR = os.path.join(HOME, "papers", "paper-howorka")
FIGDIR = os.path.join(PAPERDIR, "figures", "")

DATADIR = os.path.join(HOME, "fields")

fields.set_dir(DATADIR)


kT = physp.kT
eta = physp.eta
rMolecule = geop.rMolecule
params = dict(rMolecule=rMolecule)
eps=1e-2

l=np.linspace(rMolecule+eps,10*rMolecule,100)

def Cp(l,rMolecule):
    return 1./(1.-(9./16.)*rMolecule*(1./l))

def Cn(l,rMolecule):
    alpha=acosh(l/rMolecule)
    sum = 0.
    for n in range(1,1000):
        sum+=float(n*(n+1))/float((2*n-1))/float((2*n+3))*((2*sinh((2*n+1)*alpha)+(2*n+1)*sinh(2*alpha))/(4*(sinh((n+.5)*alpha))**2-(2*n+1)**2*(sinh(alpha))**2)-1)
    return (4./3.)*sinh(alpha)*sum
