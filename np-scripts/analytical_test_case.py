" analytical test problem to validate 2D solver "
from nanopores import *

# --- create 2D geometry ---

Rz = 2*nm
Rx = 2*nm
domain = Box([0, -Rz], [Rx, Rz])

domain.addsubdomains(fluid = domain)
domain.addboundaries(
    lowerb = domain.boundary("bottom"),
    upperb = domain.boundary("top"),
    sideb = domain.boundary("right"),
    centerb = domain.boundary("left")
)

geo = domain.create_geometry(lc=.5*nm)

# --- define physical parameters ---

phys = Physics("pore_molecule", geo)
pnps = PNPSAxisym(geo, phys)
pnps.solve()

domain.plot()

