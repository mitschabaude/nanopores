from dolfin import *
from nanopores import *
from nanopores.physics.simplepnps import *
from mysolve import hybrid_solve, newton_solve

nm = 1.

add_params(
bV = -0.05, # [V]
dnaqsdamp = .5,
h = .5,
damp = 1.,
bulkcon = 300.,
tol = 1e-15,
imax = 10,
taylorhood = False,
Rx = 8*nm,
Ry = 8*nm,
l0 = 9*nm,
iterative = False,
)
print PARAMS

geo_name = "H_geo"


geo_params = dict(
x0 = None,
boxfields = True,
Rx = Rx,
Ry = Ry,
l0 = l0
)

phys_params = dict(
Membraneqs = -0.00,
bulkcon = bulkcon,
bV = bV,
dnaqsdamp = dnaqsdamp,
)

generate_mesh(h, geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)
phys = Physics("pore", geo, **phys_params)

plot(geo.boundaries)
print "number of elements: ", geo.mesh.num_cells()
print "number of vertices: ", geo.mesh.num_vertices()

if geo.parameter("x0") is None:
    exec("from nanopores.geometries.%s.subdomains import synonymes" %geo_name)
    geo.import_synonymes({"moleculeb":set()})
    geo.import_synonymes(synonymes)
    
if not taylorhood:
    ku = 1
    beta = .01
else:
    ku = 2
    beta = .0

#print
#print "# solve pnp with fixed point method"
#pnp = PNPFixedPoint(geo, phys, cyl=True, inewton=1, ipicard=imax, verbose=True, iterative=iterative)
#t = Timer("solve")
#hybrid_solve(pnp)
#print "CPU time (solve): %s [s]" % (t.stop(),)
#pnp.visualize()

print
print "# solve pnps with fixed point method"
#SimplePoissonProblem.method["reuse"] = False
#SimpleStokesProblem.method["reuse"] = False
pnps = PNPSFixedPoint(geo, phys, cyl=True, beta=beta, ku=ku,
    inewton=1, ipicard=imax+10, tolnewton=tol, verbose=True, iterative=iterative)
t = Timer("solve")
hybrid_solve(pnps)
print "CPU time (solve): %s [s]" % (t.stop(),)
#pnps.visualize()

print "# solve with hybrid method"
SimpleStokesProblem.method["reuse"] = False
pnpsH = PNPSHybrid(geo, phys, cyl=True, beta=beta, damp=damp, ku=ku,
    inewton=1, ipicard=imax, tolnewton=tol, verbose=True, nverbose=True, iterative=iterative)
t = Timer("solve")
hybrid_solve(pnpsH)
print "CPU time (solve): %s [s]" % (t.stop(),)

print
print "# solve with newton's method"
SimplePNPSProblem.method["iterative"] = iterative
pnpsN = NonlinearPDE(geo, SimplePNPSProblem, phys=phys, cyl=True, beta=beta, ku=ku)
pnpsN.imax = imax
pnpsN.newtondamp = damp
pnpsN.tolnewton = tol
t = Timer("solve")
newton_solve(pnpsN)
print "CPU time (solve): %s [s]" % (t.stop(),)
#pnps.visualize()

v, cp, cm, u, p = pnps.solutions()
vN, cpN, cmN, uN, pN = pnpsH.solutions()
plot(v - vN)
plot(u - uN)
interactive()

# plot
pnps.estimators["err hybrid i"].name = "fixed point"
pnps.estimators["err hybrid time"].name = "fixed point"
pnpsH.estimators["err hybrid i"].name = "hybrid"
pnpsH.estimators["err hybrid time"].name = "hybrid"
pnpsN.estimators["err newton i"].name = "newton"
pnpsN.estimators["err newton time"].name = "newton"

from matplotlib import pyplot
pnps.estimators["err hybrid i"].newtonplot()
pnpsH.estimators["err hybrid i"].newtonplot(fig=False)
pnpsN.estimators["err newton i"].newtonplot(fig=False)

pnps.estimators["err hybrid time"].newtonplot()
pnpsH.estimators["err hybrid time"].newtonplot(fig=False)
pnpsN.estimators["err newton time"].newtonplot(fig=False)
pyplot.xlabel("time [s]")
pyplot.xscale("log")

#saveplots("hybrid_vs_newton", PARAMS)
showplots()



