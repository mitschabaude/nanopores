# (c) 2016 Gregor Mitscha-Baude
import dolfin, os
import numpy as np
from nanopores.models import Howorka
from nanopores.models.diffusion import friction_tensor, friction
import nanopores.tools.fields as fields
import nanopores
from nanopores import (LinearPBGoalOriented, LinearPBAxisymGoalOriented,
                       PNPSAxisym, PNPS, StokesProblem, PNPProblem, HOME,
                       Physics, user_params)
DATADIR = os.path.join(HOME, "Dropbox", "nanopores", "fields")
FIGDIR = os.path.join(HOME, "Dropbox", "nanopores", "figures")
fields.set_dir(DATADIR)

def pbpnps(geo, phys, cyl=False, frac=0.5, Nmax=1e4, cheapest=False,
                      taylorhood=False, stokesLU=False, **kwargs):
    if not stokesLU and not cyl:
        StokesProblem.method["iterative"] = True
    else:
        StokesProblem.method["iterative"] = False
    #PNPProblem.method["iterative"] = False
    PNPProblem.method["kparams"]["relative_tolerance"] = 1e-10
    PNPProblem.method["kparams"]["absolute_tolerance"] = 1e-6
    PNPProblem.method["kparams"]["nonzero_intial_guess"] = False
    #StokesProblemEqualOrder.beta = 1.
    StokesProblem.method["kparams"].update(
        monitor_convergence = False,
        relative_tolerance = 1e-10,
        absolute_tolerance = 1e-5,
        maximum_iterations = 2000,
        nonzero_initial_guess = True,
        )
    PNPS.tolnewton = 1e-3
    PNPS.alwaysstokes = True

    LinearPB = LinearPBAxisymGoalOriented if cyl else LinearPBGoalOriented
    PNPStokes = PNPSAxisym if cyl else PNPS
    z = phys.dim - 1
    bV = phys.bV
    phys.bV = 0.
    goal = (lambda v : phys.Fbare(v, z)) if geo.parameter("x0") else (lambda v : phys.CurrentPB(v))
    pb = LinearPB(geo, phys, goal=goal)
    phys.bV = bV
    pb.maxcells = Nmax
    pb.marking_fraction = frac
    if cheapest:
        pb.estimate = pb.estimate_cheap
    refined = True
    i = 0

    print "Number of cells:", pb.geo.mesh.num_cells()
    while refined:
        i += 1
        print "\nSolving PB."
        pb.single_solve()
        print "\nError estimation."
        (ind, err) = pb.estimate()
        print "\nMesh refinement."
        refined = pb.refine(ind)
        if not refined:
            print "Maximal number of cells reached."
        else:
            print "New total number of cells:", pb.geo.mesh.num_cells()

    pnps = PNPStokes(pb.geo, phys, v0=pb.solution, taylorhood=taylorhood)
    print "\nSolving PNPS."
    dofs = pnps.dofs()
    print "  Degrees of freedom: %d" % dofs
    newton_iter = pnps.newton_solve()
    print "  Newton iterations:", newton_iter
    return pb, pnps

def pnps(geo, phys, cyl=False, taylorhood=False, stokesLU=True, **kwargs):
    if not stokesLU and not cyl:
        StokesProblem.method["iterative"] = True
    else:
        StokesProblem.method["iterative"] = False
    #PNPProblem.method["iterative"] = False
    PNPProblem.method["kparams"]["relative_tolerance"] = 1e-10
    PNPProblem.method["kparams"]["absolute_tolerance"] = 1e-6
    PNPProblem.method["kparams"]["nonzero_intial_guess"] = False
    #StokesProblemEqualOrder.beta = 1.
    StokesProblem.method["kparams"].update(
        monitor_convergence = False,
        relative_tolerance = 1e-10,
        absolute_tolerance = 1e-5,
        maximum_iterations = 2000,
        nonzero_initial_guess = True,
        )
    PNPS.tolnewton = 1e-3
    PNPS.alwaysstokes = True
    PNPStokes = PNPSAxisym if cyl else PNPS

    pnps = PNPStokes(geo, phys, taylorhood=taylorhood)
    print "\nSolving PNPS."
    dofs = pnps.dofs()
    print "  Degrees of freedom: %d" % dofs
    newton_iter = pnps.newton_solve()
    print "  Newton iterations:", newton_iter
    return pnps

#@solvers.cache_forcefield("howorka_nonzero_u")
#def F(x, v=None, dim=3, **params):
#    values = []
#    setup = Howorka.setup3D if dim==3 else Howorka.setup2D
#    cyl = dim==2
#    for x0 in x:
#        geo, phys = setup(x0=x0, **params)
#        if v is not None:
#            phys.update(UMol=tuple(v))
#        #dolfin.plot(geo.submesh("solid"), key="b", title="solid mesh")
#        pb, pnps = pbpnps(geo, phys, cyl=cyl, **params)
#        #dolfin.plot(geo.submesh("solid"), key="b", title="solid mesh")
#        values.append(pnps.forces())
#        pnps.visualize("fluid")
#    F, Fel, Fdrag = tuple(zip(*values))
#    return dict(F=F, Fel=Fel, Fdrag=Fdrag)

#print F([[0.2,0.,4.6798]], Nmax=5e4, UMol=(0., 0.1, 0.1), dim=3, dnaqsdamp=0.5,
#        taylorhood=False, cheapest=True, cache=False, h3D=8.,
#        stokesLU=True)
#print F([[0.,0.,4.6798]], Nmax=5e4, UMol=(0., 0.1), dim=2, dnaqsdamp=0.5,
#        taylorhood=False, cheapest=True, cache=False, h3D=8.,
#        stokesLU=True)

class Values(list):
    def __str__(self):
        return str(self[-1])

def velocity_iteration(setup, imax=15):
    dim = setup.geop.dim
    gamma = friction_tensor(setup) if dim==3 else np.diag([friction(setup)])
    geo = setup.geo
    #dolfin.plot(geo.boundaries)
    #dolfin.interactive()
    f, v, dv = Values(), Values(), Values()

    # iteratively update v
    v.append(np.array([0.]*dim))
    for i in range(imax):
        print "\n--- Loop %d ---" %(i+1,)
        phys = Physics("howorka", geo, **setup.physp)
        phys.update(Qmol=phys.Qmol*phys.qq)
        phys.update(UMol=tuple(v[-1]))
        pde = pnps(geo, phys, cyl=setup.phys.cyl)
        Force, Fel, Fdrag = pde.forces()

        if dim==2: Force = [Force[1]]

        f.append(1e-12*np.array(Force))
        print "f =", f
        dv0 = np.linalg.solve(gamma, f[-1])
        if dim==2: dv0 = np.array([0., float(dv0)])
        dv.append(dv0)
        print "dv =", dv
        v.append(v[-1] + dv[-1])
        print "v =", v

    #pde.visualize()
    return f, v[1:], dv

@nanopores.tools.solvers.cache_forcefield("howorka_velo1")
def velocities(X, **params):
    v0, v1 = [], []
    for x0 in X:
        setup = Howorka.Setup(x0=x0, **params)
        f, v, dv = velocity_iteration(setup, 4)
        assert np.linalg.norm(1e12*f[-1]) < 1e-3
        #assert np.linalg.norm(dv[-1]) < 1e-3*np.linalg.norm(dv[0])
        v0.append(list(v[0]))
        v1.append(list(v[-1]))
    return dict(v0=v0, v1=v1)

# working 3D setup
#x = [0.2, 0., 0.]
#setup = Howorka.Setup(dim=3, Nmax=1.5e5, h=1., x0=x, dnaqsdamp=0.5, Qmol=-1.)
params = user_params(dim=2, Nmax=2e4, h=.5, dnaqsdamp=0.25,
                     Qmol=-1., bulkcon=300.)
do_plot = False

Z = np.linspace(-6., 6., 42)
X = [[0.,0.,z] for z in Z]
#X = [[0.,0.,0.]]
print velocities(X, nproc=7, **params)

if do_plot:
    x = [0., 0., 0.]

    imax = user_params(imax=15)["imax"]
    setup = Howorka.Setup(x0=x, **params)
    f, v, dv = velocity_iteration(setup, imax)

    import matplotlib.pyplot as plt
    dim = params["dim"]
    plt.semilogy(1e12*np.sqrt(np.sum(np.array(f)**2, 1)), "s-")
    plt.ylabel("force on molecule [pN]")
    plt.xlabel("# iterations")
    #plt.figure()
    #plt.plot(np.array(v)[1:,dim-1], "s-")
    #plt.ylabel("infered velocity [m/s]")
    #plt.xlabel("# iterations")
    nanopores.savefigs("howorka_velocity_z4.6", FIGDIR)