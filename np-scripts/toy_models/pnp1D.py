" 1D PNP "
from dolfin import *
from nanopores import *
from nanopores.physics.simplepnps import *

# --- create 1D geometry ---

h = 20.
hmem = 2.

domain = Box([-h/2], [h/2])
membrane = Box([-hmem/2], [hmem/2])
lowerb = domain.boundary("left")
upperb = domain.boundary("right")

domain.addsubdomains(
    fluid = domain - membrane,
    membrane = membrane
)
domain.addboundaries(
    lowerb = lowerb,
    upperb = upperb,
    chargedmembraneb = membrane.boundary(),
)
domain.params = dict(
    lscale = 1e9,
)
domain.synonymes = dict(
    solid = "membrane",
)

geo = domain.create_geometry(lc=.01)  
print geo

# --- define physical parameters ---

phys_params = dict(
Membraneqs = -0.0,
bulkcon = 1e2,
bV = -.1,
dnaqsdamp = .25
)
phys = Physics("pore", geo, **phys_params)
#print phys.default_synonymes
#geo.import_synonymes(phys.default_synonymes)

print geo

# --- define and solve PDE ---

pnps = NonlinearPDE(geo, SimplePNPProblem, phys=phys, axisymmetric=False)
pnps.imax = 20
pnps.newtondamp = 1.
pnps.maxcells = 5e4
t = Timer("solve")
pnps.solve(refinement=False)
print "CPU time (solve): %s [s]" % (t.stop(),)
pnps.visualize()
#domain.plot()

