from dolfin import *
from nanopores import *
from nanopores.physics.simplepnps import *

geo_name = "H_geo"
nm = 1e-9

geo_params = dict(
x0 = None,
boxfields = True,
#Rx = 300*nm,
#Ry = 30*nm,
)

phys_params = dict(
Membraneqs = -0.01,
bulkcon = 3e2,
bV = -.1,
dnaqsdamp = .25
)

generate_mesh(.5, geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)
phys = Physics("pore", geo, **phys_params)

plot(geo.subdomains)
plot(geo.boundaries)
print geo

if geo.parameter("x0") is None:
    exec("from nanopores.geometries.%s.subdomains import synonymes" %geo_name)
    geo.import_synonymes({"moleculeb":set()})
    geo.import_synonymes(synonymes)
"""
# solve with hybrid method
pnps = NonlinearPDE(geo, SimplePNPProblem, phys=phys, cyl=True)
pnps.imax = 20
pnps.newtondamp = 1.
t = Timer("solve")
pnps.solve()
print "CPU time (solve): %s [s]" % (t.stop(),)
print phys
pnps.visualize()
"""
# solve with newton's method
pnps = NonlinearPDE(geo, SimplePNPSProblem, phys=phys, cyl=True)
pnps.imax = 20
pnps.newtondamp = 1.
t = Timer("solve")
pnps.solve()
print "CPU time (solve): %s [s]" % (t.stop(),)
pnps.visualize()
"""
tol = 1e-2
damp = 1.
S = pnps.solvers.values()[0] 
S.newtondamp = damp

for i in range(20):
    #plot(pnps.functions.values()[0].sub(0)) # for debugging
    #interactive()
    pnps.visualize()
    S.solve()
    print 'Relative L2 Newton error:',S.relerror()
    if S.convergence(tol):
        print 'linf Norm of Newton update:', \
                    norm(S.problem.u.vector(),'linf'), \
                    '<=', tol ,' \n  ==> break loop \n'
        break
print "Newton iterations:",i+1
print 'Relative L2 Newton error:',S.relerror()

pnps.visualize()

"""
pnps.print_results()
#pnps.estimators["est"].plot(rate=-.5)

