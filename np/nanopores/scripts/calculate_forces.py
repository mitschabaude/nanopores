''' calculate forces on molecule depending on midpoint '''

from nanopores import *
from dolfin import *
import sys, argparse, math

# general parameters
rMolecule = 0.55 # [nm]
Qmol = -2. # [q*C]
dnaqsdamp = 0.1
bV0 = 0.01
z0 = 7.5 # [nm]


# geo parameters 3D
geo_name = "H_cyl_geo"
nm = import_vars("nanopores.geometries.%s.params_geo" %geo_name)["nm"]
params = dict(
rMolecule = 0.55*nm,
#moleculeblayer = True,
)

# geo parameters 2D
geo_name2D = "H_geo"
nm = 1e-9
params2D = dict(
rMolecule = 0.55*nm,
moleculeblayer = True,
boxfields = True,
)

# physical parameters (2D + 3D)
phys_params = dict(
Membraneqs = -0.,
#bV = 0.01,
Qmol = Qmol*qq,
bulkcon = 3e2,
#uppermbias = 1.,
#lowermbias = -.02,
dnaqsdamp = dnaqsdamp,
couplebVtoQmol = True,
bV0 = bV0,
)


IllposedNonlinearSolver.newtondamp = 1.
#set_log_level(PROGRESS)

# solver parameters 3D
default_maxcells = 16e4
PNPS.imax = 50
PNPS.tolnewton = 1e0
#PNPProblem.method["kparams"]["monitor_convergence"] = True
PNPProblem.method["kparams"]["maximum_iterations"] = 600
#PNPProblem.method["kparams"]["error_on_nonconvergence"] = False
#PNPProblem.method["iterative"] = False
PNPProblem.method["kparams"]["absolute_tolerance"] = PNPS.tolnewton*1e-3
PNPProblem.method["kparams"]["relative_tolerance"] = 1e-8
StokesProblemEqualOrder.method["iterative"] = False
StokesProblemEqualOrder.method["kparams"]["absolute_tolerance"] = PNPS.tolnewton*1e-3
#StokesProblemEqualOrder.method["kparams"]["monitor_convergence"] = True
#LinearPBProblem.method["kparams"]["monitor_convergence"] = True
PoissonProblem.method["kparams"]["relative_tolerance"] = 1e-6
LinearPBProblem.method["kparams"]["relative_tolerance"] = 1e-6

# solver parameters 2D
default_maxcells2D = 12e4
PNPSAxisym.tolnewton = 1e0


def calculate_forces(x0, pid="", clscale=10.0, refinement=True, maxcells=default_maxcells):
    ''' calculate forces on molecule depending on midpoint '''
    
    t = Timer("Mesh Generation")
    generate_mesh(clscale, geo_name, pid=pid, x0=x0, **params)
    meshfile = "/".join([DATADIR, geo_name, "mesh", "mesh%s.xml" %pid])
    geo = geo_from_name(geo_name, mesh=Mesh(meshfile), x0=x0, **params)
    phys = Physics("pore_molecule", geo, **phys_params)
    print "CPU Time (mesh generation):",t.stop()
    print "hmin:", geo.mesh.hmin()
    
    t = Timer('PB')
    goal = lambda v : phys.Fbare(v, 2) + phys.CurrentPB(v)
    pb = LinearPBGoalOriented(geo, phys, goal=goal)
    pb.maxcells = maxcells
    pb.marking_fraction = 0.2
    pb.solve(refinement=refinement)
    print "CPU Time (PB):",t.stop()

    t = Timer('PNPS')
    geo = pb.geo
    v0 = pb.solution
    #pb.visualize("solid")
    pnps = PNPS(geo, phys, v0=v0)
    #pnps = PNPS(geo, phys)
    i = pnps.solve(visualize=False)
    while i==50:
        print "\nRestarting Newton iteration!"
        pnps.__init__(geo, phys)
        pnps.solvers["PNP"].newtondamp *= 0.8
        i = pnps.solve()
    print "Newton iterations:",i

    print "CPU Time (PNPS):",t.stop()
    pnps.print_results()
    f = pnps.get_functionals()
    if any(math.isnan(f[s]) for s in f):
        raise Exception("NaN occured in force dict, probably caused by PETSc failure.")
    return pnps.get_functionals()
    
    
def calculate_forces2D(x0, pid="", clscale=.8, refinement=True, maxcells=default_maxcells2D):
    ''' calculate forces on molecule depending on midpoint '''
    nm = 1e-9 # by convention, nm == 1. in mesh generation script
    x0 = map(lambda x:x*nm, x0)
    
    t = Timer("Mesh Generation")
    generate_mesh(clscale, geo_name2D, pid=pid, x0=x0, **params2D)
    meshfile = "/".join([DATADIR, geo_name2D, "mesh", "mesh%s.xml" %pid])
    geo = geo_from_name(geo_name2D, mesh=Mesh(meshfile), x0=x0, **params2D)
    phys = Physics("pore_molecule", geo, **phys_params)
    print "CPU Time (mesh generation):",t.stop()
    print "hmin:", geo.mesh.hmin()
    
    t = Timer('PB')
    goal = lambda v : phys.Fbare(v, 1) + phys.CurrentPB(v)
    pb = LinearPBAxisymGoalOriented(geo, phys, goal=goal)
    pb.maxcells = maxcells
    pb.marking_fraction = 0.5
    pb.solve(refinement=refinement)
    print "CPU Time (PB):",t.stop()

    t = Timer('PNPS')
    geo = pb.geo
    v0 = pb.solution
    pnps = PNPSAxisym(geo, phys, v0=v0)
    #pnps = PNPSAxisym(geo, phys)
    i = pnps.solve(visualize=False)
    while i==50:
        print "\nRestarting Newton iteration!"
        pnps.__init__(geo, phys)
        pnps.solvers["PNP"].newtondamp *= 0.8
        i = pnps.solve()
    print "Newton iterations:",i

    print "CPU Time (PNPS):",t.stop()
    pnps.print_results()
    
    # make forces 3D
    f = pnps.get_functionals()
    for s in {"Fp", "Fshear", "Fbare", "Fbarevol"}:
        f[s+"2"] = f[s+"1"]
        f[s+"1"] = 0.
        
    if any(math.isnan(f[s]) for s in f):
        raise Exception("NaN occured in force dict, probably caused by PETSc failure.")
    return f
    
def calculate2D(clscale=.8, refinement=True, maxcells=10e4, pid="", **params):
    globals().update(params)
    nm = 1e-9
    global params2D, phys_params
    params2D["rMolecule"] = rMolecule*nm
    if "x0" in params:
        x0 = params["x0"]
    else:
        x0 = [0.,0.,z0*nm]
    phys_params.update(dict(
        Qmol = Qmol*qq,
        dnaqsdamp = dnaqsdamp,
        bV0 = bV0,
    ))
    forces = calculate_forces2D(x0, pid=pid, clscale=clscale, refinement=refinement, maxcells=maxcells)
    F = sum(forces[s] for s in ["Fbarevol2", "Fp2", "Fshear2"])
    return F
    
    
if __name__ == "__main__":
    # parse user arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('x0', default="[0.0, 0.0, 0.0]", help='Molecule position')
    parser.add_argument('pid', default="", help='Process ID for .out .msh .xml files')
    parser.add_argument('clscale', default=15.0, type=float, help='Scale')
    parser.add_argument('dim', default=3, type=int, help='Dimension')
    args, unknown = parser.parse_known_args()
    #print eval(args.x0), args.pid, args.clscale
    print
    print "dim =",args.dim
    if args.dim==2:
        print calculate_forces2D(eval(args.x0), args.pid, args.clscale)
    elif args.dim==3:
        print calculate_forces(eval(args.x0), args.pid, args.clscale)

