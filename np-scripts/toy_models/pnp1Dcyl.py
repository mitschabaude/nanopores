" 1D cylindrical PB, modelling crossection of pore "
# Comment: this is a prime example of simplicity and flexibility :)
# the file solves a real PDE with precise specifications but depends ONLY on core library functions!

from dolfin import *
from nanopores import *
from nanopores.physics.simplepnps import *

# --- create 1D geometry ---

R = 1. # [nm] pore radius

domain = Interval(0., R)
domain.addsubdomain(domain, "fluid")
domain.addboundaries(
    wall = domain.boundary("right"),
    center = domain.boundary("left")
)
domain.params["lscale"] = 1e9
domain.synonymes = dict(
    water = "fluid",
)

geo = domain.create_geometry(lc=.01)  
print geo

# --- define physical parameters ---

phys_params = dict(
    surfcharge = dict(wall = -0.1*qq/nm**2),
    bulkcon = 300,
    #volcharge = -.01*qq/nm**3, # to account for \partial_{zz}\phi?
)
phys = Physics("electrolyte", geo, **phys_params)
print geo

# --- solve pdes ---

# nonlinear PB
pb = solve_pde(SimplePBProblem, geo, phys, cyl=True, iterative=False, tolnewton=1e-10)
# linear PB -- imax=1 implies only one Newton iteration
pblin = solve_pde(SimplePBProblem, geo, phys, cyl=True, iterative=False, imax=1)

u0 = pb.solution
u1 = pblin.solution
plot1D({"linear PB":u1, "nonlinear PB":u0},
        (0., R, 101), "x", dim=1, axlabels=("r [nm]", "potential [V]"))
showplots()

# PNP doesn't converge and probably doesn't make sense, since bulk concentration of ions can not be assigned
#pnps = solve_pde(SimplePNPProblem, geo, phys, cyl=True, visualize=True)

