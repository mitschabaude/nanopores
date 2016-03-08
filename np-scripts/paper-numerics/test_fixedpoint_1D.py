from dolfin import *
from nanopores import *
from nanopores.physics.simplepnps import *
from mysolve import hybrid_solve, newton_solve

add_params(
bV = -0.1, # [V]
rho = -0.05,
bulkcon = 300.,
imax = 10,
)

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
geo = domain.create_geometry(lc=.1)

phys_params = dict(
    Membraneqs = rho,
    bulkcon = bulkcon,
    v0 = dict(upperb = 0., lowerb = bV),
)
phys = Physics("pore", geo, **phys_params)

# --- define and solve PDE ---
#p = solve_pde(SimplePoissonProblem, geo, phys)

pnp = PNPFixedPoint(geo, phys, inewton=1, ipicard=imax, verbose=True)
#t = Timer("solve")
hybrid_solve(pnp)
#print "CPU time (solve): %s [s]" % (t.stop(),)
#pnp.visualize()

v, cp, cm = pnp.solutions()
#cp, cm, v = pnp.solutions()
#v  = p.solution
plot1D({"potential": v}, (-h/2, h/2, 101), "x", dim=1, axlabels=("z [nm]", "potential [V]"))
plot1D({"c+": cp, "c-":cm},  (hmem/2, h/2, 101), "x", dim=1, axlabels=("z [nm]", "concentrations [mol/m^3]"))       
showplots()

