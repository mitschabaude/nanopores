from nanopores import *
from dolfin import *

geo_name = "P_geo"
geo_params = {"x0": [0.,0.,0.]}  # any parameter change here will not affect mesh generation yet
phys_params = {"bV":0.1, "MembraneSiNqs":-30e-3, "bulkcon":1e3, "DNAqs": -0.1362}

clscale = 1.
meshgen_dict = generate_mesh(clscale, geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)
phys = Physics("pore_dna", geo, **phys_params)
print phys.charge

PNPSAxisym.newtondamp=1
PNPSAxisym.imax=20

phys.bV = 0.0
#goal = lambda v : phys.Fbare(v, 1) + phys.CurrentPB(v)
#pb = LinearPBAxisymGoalOriented(geo, phys, goal=goal)
pb = LinearPBAxisym(geo, phys)
pb.maxcells = 4e4
pb.marking_fraction = 0.5
pb.solve(refinement=True)
pb.estimate()
#pb.estimators["err"].plot(rate=-1./2.)
#pb.visualize()

geo = pb.geo
phys.bV = phys_params["bV"]

pnps = PNPSAxisym(geo, phys)
pnps.maxcells = 2e4
pnps.marking_fraction = 0.5
#pnps.solvers.pop("Stokes")
pnps.StokesProblem = StokesProblemAxisym
pnps.solve(refinement=False, print_functionals=True)
#pnps.print_results()

(v,cp,cm,u,p) = pnps.solutions(deepcopy=True)
print "phys.Fbare [pN]:",1e12*assemble(phys.Fbare(v, 1))
pnps.visualize("fluid")

adddict = dict(
    hmin = pnps.geo.mesh.hmin(),
    N = pnps.geo.mesh.num_cells(),
    clscale = clscale
)
qoi = [adddict["N"], adddict["hmin"], pnps.get_functionals()["Feff"]]

print phys_params
print qoi
#interactive()

l0 = 50*nm
r0 = 3*nm
I = pnps.get_functionals()["Javgctr"]
V = v([r0, -l0]) - v([r0, l0])
print "I (current through pore center):",I,"[pA]"
print "V (transmembrane potential):",V,"[V]"
print "conductance I/V:",I/V,"[pS]"
print
