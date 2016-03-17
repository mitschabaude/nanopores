" 1D PNP, modelling reservoirs and membrane far away from pore "

from dolfin import *
from nanopores import *
from nanopores.physics.simplepnps import *

# --- create 1D geometry ---

h = 20.
hmem = 3.

domain = Interval(-h/2, h/2)
membrane = Interval(-hmem/2, hmem/2)
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
domain.params["lscale"] = 1e9
domain.synonymes = dict(
    solid = "membrane",
    bulkfluid = "fluid",
    pore = set()
)
geo = domain.create_geometry(lc=.01)  
print geo

# --- define physical parameters ---

phys_params = dict(
    Membraneqs = -0.0,
    bulkcon = 300.,
    bV = -.1,
)
phys = Physics("pore", geo, **phys_params)
print geo # this tests how the default synonymes in physics/pore are incorporated

# --- define and solve PDE ---

pnps = solve_pde(SimplePNPProblem, geo, phys)
v, cp, cm = pnps.solutions()

plot1D({"potential": v}, (-5., 5., 101), "x", dim=1, axlabels=("z [nm]", "potential [V]"))
plot1D({"c+": cp, "c-":cm},  (hmem/2, h/2, 101), "x", dim=1, axlabels=("z [nm]", "concentrations [mol/m^3]"))       
showplots()

