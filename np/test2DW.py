from nanopores import *
from dolfin import *

geo_name = "W_2D_geo"
geo_dict = import_vars("nanopores.geometries.%s.params_geo" %geo_name)
physical_dict = import_vars("nanopores.physics.params_physical")
default_dict = dict(geo_dict = geo_dict, physical_dict = physical_dict)

nm = geo_dict["nm"]
lsam = geo_dict["lsam"]
geo_params = {"x0": None, #[0,0,0],
              "outerfrac":0.2,
              "lsam":lsam, "r0":13*nm-lsam, "angle":40}
phys_params = {"SiNqs":-0.022, "SAMqs":-0.078, "bV": 0.2}

meshgen_dict = generate_mesh(0.5, geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)
#mesh = Mesh("/".join([DATADIR, geo_name, "mesh", "last_adapted_mesh.xml"]))
#geo = geo_from_name(geo_name, mesh=mesh, **params)
print "Number of cells:",geo.mesh.num_cells()
print "hmin:", geo.mesh.hmin()

plot(geo.subdomains, title="geometry")

PNPSAxisym.imax = 50
PNPSAxisym.maxcells = 2e4
PNPSAxisym.marking_fraction = 0.5
PNPSAxisym.tolnewton = 1e-2

pnps = PNPSAxisym(geo, **phys_params)

pnps.solve(refinement=True, save_mesh=False, visualize=False)
#pnps.solvers.pop("Stokes")
pnps.print_results()
plot(geo.boundaries, title="final boundaries")


'''
for est in pnps.estimators.values():
    print "\n%s estimator:\n" %est.name, est[:]
    est.save_to_matlab()

#print "Convergence rates:\n",pnps.estimators["h1"].rates()

#(v,cp,cm,u,p) = pnps.solutions()
#submesh = pnps.geo.submesh("porecenter")
#for f in {v,cp,cm,u}:
#    adaptfunction(f, submesh, assign=True)
#Jp = (-D*grad(cp) - mu*cp*grad(v) + cp*u)[1]
#plot(Jp)
'''
pnps.visualize("fluid")
