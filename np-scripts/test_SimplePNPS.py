from dolfin import *
from nanopores import *
from nanopores.physics.simplepnps import *

geo_name = "H_geo"
nm = 1e-9
z0 = 7.5*nm

geo_params = dict(
x0 = None,
#x0 = [0.,0.,z0],
#x0 = [0., 0., -8.372*nm],
#rMolecule = 0.4*nm,
#lcMolecule = nm*0.1,
#moleculeblayer = True,
boxfields = True,
#Rx = 300*nm,
#Ry = 30*nm,
)

phys_params = dict(
Membraneqs = -0.0,
bulkcon = 3e2,
bV = -.9,
dnaqsdamp = .5
)

generate_mesh(.9, geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)
phys = Physics("pore", geo, **phys_params)

plot(geo.subdomains)
plot(geo.boundaries)

if geo.parameter("x0") is None:
    exec("from nanopores.geometries.%s.subdomains import synonymes" %geo_name)
    geo.import_synonymes({"moleculeb":set()})
    geo.import_synonymes(synonymes)

pnps = NonlinearPDE(geo, SimplePNPProblem, phys=phys, axisymmetric=True)
pnps.imax = 20
pnps.newtondamp = 1.

tol = 1e-2
damp = 1.
S = pnps.solvers.values()[0] 
S.newtondamp = damp

for i in range(20):
    #plot(pnps.functions.values()[0].sub(0)) # for debugging
    #interactive()
    S.solve()
    print 'Relative L2 Newton error:',S.relerror()
    if S.convergence(tol):
        print 'linf Norm of Newton update:', \
                    norm(S.problem.u.vector(),'linf'), \
                    '<=', tol ,' \n  ==> break loop \n'
        break
print "Newton iterations:",i+1
print 'Relative L2 Newton error:',S.relerror()

pnps.print_results()
pnps.visualize()

