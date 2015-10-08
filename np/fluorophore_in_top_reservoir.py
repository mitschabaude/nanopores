''' this script is designed to determine the electric field induced near the nanopore by a cylindrical upper liquid reservoir filled with charged particles (fluorophores), given the concentration, single molecule charge, and dimensions of the reservoir.
the charge is smeared out across the liquid.
'''

from nanopores import *
from dolfin import *

# ----- physical parameters -----
mm = 1e-3
lx = 10*mm # radius [m]
ly = 20*mm # height [m]
uM = 1e-3
bulkcon = 10*uM # bulk concentration [mM] [mol/m**3]
Qmol = -3*qq # charge per molecule

# ----- discretization parameters -----
N = 40 # number of points in every space direction

# volume charge [C/m**3]
volcharge = bulkcon*mol*Qmol
print "volume charge [C/m**3]:",volcharge

mesh = RectangleMesh(Point(0.,0.), Point(lx, ly), N, N)

geo = Geometry(
    mesh = mesh,
    synonymes = {"fluid":"all"})
    
phys = Physics(
    geo = geo,
    surfcharge = {},
    volcharge = volcharge,)
    
PoissonProblemPureNeumannAxisym.method["iterative"] = False
pde = LinearPDE(geo, PoissonProblemPureNeumannAxisym)
pde.solve()
pde.visualize()

