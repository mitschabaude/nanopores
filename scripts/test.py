from nanopores import *
from dolfin import *

geo_name = "H_cyl_geo"
nm = import_vars("nanopores.geometries.%s.params_geo" %geo_name)["nm"]
z0 = 4.0*nm

geo_params = dict(
x0 = [0.,0.,z0],
R = 15*nm,
Rz = 15*nm,
rMolecule = 0.4*nm,
r0 = 1.*nm,
moleculeblayer = True,
)

phys_params = {"Membraneqs": -0.03, "bV": 0.001, "Qmol":-3.*qq, "bulkcon":1000}

clscale = 8.0
mesh_dict = generate_mesh(clscale, geo_name, **geo_params)
geo = geo_from_name(geo_name, **geo_params)
phys = Physics("pore_molecule", geo, **phys_params)

# subd = "fluid"
# submesh = geo.submesh(subd)
# geo_sub = geo_from_subdomains(submesh, "nanopores.%s.subdomains" %geo.params["name"], **geo.params)
# plot(geo_sub.boundaries, title=("boundaries on %s" %subd), elevate=-3e1)
# plot(submesh, title=("initial mesh on %s" %subd), wireframe=True, elevate=-3e1)
# interactive()
# exit()

IllposedLinearSolver.stab = 1.0/nm
PNPProblem.dnaqsdamp = 0.7
PNPProblem.method["kparams"]["maximum_iterations"] = 10000
#PNPProblem.method["kparams"]["iterative"] = False
#StokesProblemEqualOrder.method["kparams"]["monitor_convergence"] = True

PNPS.imax = 50
PNPS.tolnewton = 1e0
StokesProblemEqualOrder.method["kparams"]["absolute_tolerance"] = PNPS.tolnewton*1e-3

#set_log_level(PROGRESS)

#IllposedNonlinearSolver.newtondamp = 0.8
print("hmin:", geo.mesh.hmin())
init_cells = geo.mesh.num_cells()
print("number of cells:", init_cells)

timer = Timer('Start')
pnps = PNPS(geo, phys)
stokes = pnps.solvers.pop("Stokes")
pnps.maxcells = 5*init_cells
pnps.marking_fraction = 0.2
pnps.solve(refinement=True, save_mesh=False, visualize=False, print_functionals=False)

pnps.print_results()
print("All Functionals are multiplied by", PNPS.Functional_mult)

print("Total time: ", timer.stop())

Jl = [pnps.get_functionals()[j] for j in ["Javgtop","Javgctr","Javgbtm"]]
qdict = dict(
    hmin = geo.mesh.hmin(),
    N = geo.mesh.num_cells(),
    clscale = clscale,
    solvers = list(pnps.solvers.keys()),
    J = max(Jl),
)
#qoi = [dqict[i] for i in gdict.items()
print(list(qdict.items()))
for est in list(pnps.estimators.values()):
    est.plot(rate=-0.3333)

#list_timings()
pnps.visualize("pore")
