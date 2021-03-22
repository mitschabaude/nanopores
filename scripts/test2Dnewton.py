from nanopores import *
from dolfin import *

geo_name = "H_geo"
geo_params = {"x0":None} #import_vars("params_nice_2D")
phys_params = {"bV": 0.1, "bulkcon":300, "dnaqsdamp":0.1}

meshgen_dict = generate_mesh(1., geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)
#mesh = Mesh("/".join([DATADIR, geo_name, "mesh", "last_adapted_mesh.xml"]))
#geo = geo_from_name(geo_name, mesh=mesh, **geo_params)
print("Number of cells:",geo.mesh.num_cells())
# plot(geo.indicator("solid"), title="geometry")
# plot(geo.mesh, title="initial mesh")

if geo.parameter("x0") is None:
    exec("from nanopores.geometries.%s.subdomains import synonymes" %geo_name)
    geo.import_synonymes({"moleculeb":set()})
    geo.import_synonymes(synonymes)

PNPS.imax = 16
PNPS.maxcells = 100e3
PNPS.marking_fraction = 0.5
PNPS.tolnewton = 1e-16
PNPS.alwaysstokes = True
StokesProblemAxisym.method["iterative"] = False
#PNPSProblemAxisym.method["iterative"] = True

phys = Physics("pore_molecule", geo, **phys_params)

pnps = PNPSAxisymNewton(geo, phys)
#pnps = PNPSAxisym(geo, phys)

pnps.solve(refinement=False, save_mesh=False, visualize=False)
print(phys)
pnps.print_results()


for est in list(pnps.estimators.values()):
    print("\n%s estimator:\n" %est.name, est[:])
    est.newtonplot()

# print "Convergence rates:\n",pnps.estimators["h1"].rates()

# (v,cp,cm,u,p) = pnps.solutions()
# submesh = pnps.geo.submesh("porecenter")
# for f in {v,cp,cm,u}:
#    adaptfunction(f, submesh, assign=True)
# Jp = (-D*grad(cp) - mu*cp*grad(v) + cp*u)[1]
# plot(Jp)
# list_timings()
pnps.visualize()

import resource
reskB = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print("resource consumption [MB]: ", reskB/1024.)

showplots()
