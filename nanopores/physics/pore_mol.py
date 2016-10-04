# (c) 2016 Gregor Mitscha-Baude
"""parameters for pore with spherical molecule(s) inside
(3D or axisymmetric 2D geometries are assumed)"""

from nanopores.physics.pore import *

Qmol = 0. # molecule charge [q]
Qmolq = lambda Qmol, qq: Qmol*qq
qTarget = Qmolq
rpermMol = rpermDNA
permMol = lambda rpermMol: eperm*rpermMol

def r2pi(dim):
    return dolfin.Expression("2*pi*x[0]") if dim==2 else dolfin.Constant(1.)

volcharge.update(
    molecule = "Moleculeqv",
)
permittivity.update(
    molecule = "permMol",
)

# geometry-dependent parameters
def Moleculeqv(geo, Qmolq, lscale, r2pi): # 
    "Molecule volume charge density [C/nm**3], adapted to discrete volume"
    scale = dolfin.Constant(1.0/lscale**3)
    r = scale*r2pi
    def compute(geo):
        vol = dolfin.assemble(r*geo.dx("molecule"))
        return Qmolq/vol if vol > 0. else 0.
    const = geo.constant("Moleculeqv", compute)
    return const
    
def rTarget(geo, lscale):
    return geo.params["rMolecule"]/lscale

# functionals of various continuous quantities
def Fbare(geo, r2pi, Moleculeqv, grad):
    def Fbarevol(v, i):
        dx = geo.dx("molecule")
        #scale = dolfin.Constant(lscale**(-3))
        return Moleculeqv*(-r2pi*grad(v)[i])*dx
    return Fbarevol
    
def CurrentPB(geo, r2pi, bulkcon, mu, rDPore, UT, lscale, cFarad):
    "approximation of current at 100 mV as linear functional of PB solution"
    bV0 = 0.1
    def J0(v):
        L = geo.params["lpore"]/lscale
        E = bV0/L
        Jz = 2*cFarad*bulkcon*mu*rDPore*v/UT*E* r2pi/L*geo.dx("pore")
        return Jz
    return J0
    