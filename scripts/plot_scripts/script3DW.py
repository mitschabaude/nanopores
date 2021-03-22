from nanopores import *
from dolfin import *
import json

# MPI version, run with: mpirun -n processes python _.py
parameters["std_out_all_processes"] = False;
comm = mpi_comm_world()
mpiSize = MPI.size(comm)

geo_name = "W_3D_geo"
geo_dict = import_vars("nanopores.%s.params_geo" %geo_name)
physical_dict = import_vars("nanopores.physics.params_physical")
default_dict = dict(geo_dict = geo_dict, physical_dict = physical_dict)
print("Default parameters: \n", json.dumps(default_dict, indent=4, sort_keys=True))
nm = geo_dict["nm"]

# no sam and no molecule causes a "found no facets..."
lsam = geo_dict["lsam"]
#config Fig 1f: sam, dp=26nm, bV=0.2V, c0 = 1M, pH=7.6  -> J~25nA
config1f = [{"x0": None, "lsam": lsam, "r0": 13*nm -lsam,},
            {"bV": 0.18, "bulkcon":1e3, "SAMqs": -0.031, "SiNqs": -0.022},]
#config Fig SI2b: dp=23nm, c0 = 400mM, pH=9
configs2b= [{"x0": None, "lsam": lsam, "r0": 11.5*nm -lsam,},
            {"bulkcon":4e2, "SAMqs": -0.078, "SiNqs": -0.022},]

configl = configs2b
configl.append({"membraneblayer": False,"angle":20.0})
defparams = dict((k, v) for d in configl for k, v in list(d.items()))

if mpiSize == 1:
    tmesh = Timer('Start')
    mesh_dict = generate_mesh(5.0, geo_name, **defparams)
    mesh = None #Mesh("%s/%s/mesh/mesh_bl_138k.xml" %(DATADIR,geo_name))
    tmesh.stop()
    print("Mesh generation time: ", tmesh.value())
geo = geo_from_name(geo_name, mesh=mesh, **defparams)
plot(geo.boundaries, interactive=True)

IllposedLinearSolver.stab = 1.0/nm
PNPProblem.method["iterative"] = True
PNPProblem.method["kp"] = "hypre_euclid"
PNPProblem.method["kparams"]["relative_tolerance"] = 1e-2
PNPProblem.method["kparams"]["maximum_iterations"] = 1000
StokesProblem.method["iterative"] = True

print("hmin:", geo.mesh.hmin())
print("number of cells:", geo.mesh.num_cells())

PNPS.tolnewton = 1e-2
PNPS.imax = 50

t = Timer('Start')
pnps = PNPS(geo,**defparams)
pnps.solvers.pop("Stokes")
pnps.solve(refinement=False, save_mesh=False, visualize=False, print_functionals=False)
pnps.print_results()
print("All Functionals are multiplied by", PNPS.Functional_mult)

t.stop()
print("Time for solving: ", t.value())

debye_length = 1.0/sqrt(2*physical_dict["cFarad"]*defparams["bulkcon"]/\
                        (physical_dict["rpermw"]*physical_dict["eperm"]*physical_dict["UT"]))
print("Debye length: ", debye_length)
print(defparams)

Jl = [pnps.get_functionals()[j] for j in ["Javgtop","Javgctr","Javgbtm"]]
qdict = dict(
    hmin = geo.mesh.hmin(),
    N = geo.mesh.num_cells(),
    solvers = list(pnps.solvers.keys()),
    Jl = Jl,
    Jmax = max(Jl),
)
print(qdict)

'''
for est in pnps.estimators.values():
    print "\n%s estimator:\n" %est.name, est[:]
    est.save_to_matlab()

# plot velocitiy *direction*
parameters["form_compiler"]["quadrature_degree"] = 5
u = pnps.solutions("Stokes")[0]
mesh = geo.submesh("pore")
u = adaptfunction(u, mesh)
V = StokesProblem.space(mesh).sub(0).collapse()
u0 = conditional(sqrt(u**2) > 1e-4, u/sqrt(u**2), Constant((0.,0.,0.)))
plot(u0, title="velocity direction")
'''
pnps.visualize()
