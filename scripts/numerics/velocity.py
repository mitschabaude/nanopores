# (c) 2016 Gregor Mitscha-Baude
import dolfin, os
import numpy as np
import matplotlib.pyplot as plt
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
    def tolist(self):
        return [x.tolist() for x in self]

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
    v.pop(0)
    return f, v, dv

@nanopores.tools.solvers.cache_forcefield("howorka_velo2")
def velocities(X, **params):
    v0, v1 = [], []
    for x0 in X:
        setup = Howorka.Setup(x0=x0, **params)
        f, v, dv = velocity_iteration(setup, 6)
        assert np.linalg.norm(1e12*f[-1]) < 1e-3
        #assert np.linalg.norm(dv[-1]) < 1e-3*np.linalg.norm(dv[0])
        v0.append(list(v[0]))
        v1.append(list(v[-1]))
    return dict(v0=v0, v1=v1)

def velo2force_2D(v, setup):
    geo = setup.geo
    phys = Physics("howorka", geo, **setup.physp)
    phys.update(Qmol=phys.Qmol*phys.qq)
    vv = [0.]*setup.geop.dim
    vv[-1] = v
    phys.update(UMol=tuple(vv))
    pde = pnps(geo, phys, cyl=setup.phys.cyl)
    v, cp, cm, u, p = pde.solutions()
    dolfin.plot(u)
    Force, Fel, Fdrag = pde.forces()
    return 1e-12*Force[1]

def velo2force_3D(v, setup):
    geo = setup.geo
    phys = Physics("howorka", geo, **setup.physp)
    phys.update(Qmol=phys.Qmol*phys.qq)
    phys.update(UMol=tuple(v))
    pde = pnps(geo, phys, cyl=setup.phys.cyl)
    #v, cp, cm, u, p = pde.solutions()
    #dolfin.plot(u)
    pde.visualize("fluid")
    Force, Fel, Fdrag = pde.forces()
    return 1e-12*Force[-1] #[1e-12*f for f in Force]

#@nanopores.tools.solvers.cache_forcefield("howorka_veloforce")
def nonzero_velocities_2D(V, **params):
    setup = Howorka.Setup(**params)
    gamma = friction(setup)
    print "friction gamma", gamma
    # determine F(0), only once
    if not 0. in V:
        V.append(0.)
        V.sort()

    F = [None]*len(V)
    i0 = V.index(0.)
    F[i0] = velo2force_2D(0., setup)
    F0 = F[i0]
    print "F(0)", F0

    for i, v in enumerate(V):
        if not i == i0:
            print "\n--- Velocity %d ---" %(i+1,)
            F[i] = velo2force_2D(v, setup)
            print "Velocity", v
            print "Force (exact)", F[i]
            print "Force (linear)", F0 - gamma*v

    return F, gamma, F0

params = user_params(dim=3, Nmax=1.5e5, h=1., dnaqsdamp=0.25,
                     x0=[0.2,0.,4.01], Qmol=-1., bulkcon=300.)
setup = Howorka.Setup(**params)
setup.prerefine()
velo2force_3D([0., 0.1, 0.2], setup)

do_v2f = False
redo_v2f = False
if do_v2f:
    if redo_v2f:
        params = user_params(dim=2, Nmax=2e4, h=.5, dnaqsdamp=0.25,
                             x0=[0.,0.,4.5], Qmol=-1., bulkcon=300.)
        V = list(np.linspace(-1., 1., 3))
        F, gamma, F0 = nonzero_velocities_2D(V, **params)
        fields.save_entries("howorka_velo2force_3", params, V=V, F=F, gamma=gamma, F0=F0)
        fields.update()

    dolfin.interactive()
    data = fields.load_file("howorka_velo2force_3")
    V, F, F0, gamma = tuple(data[key] for key in ["V", "F", "F0", "gamma"])
    ax = plt.axes()
    ax.axhline(y=0, color='#999999', linewidth=0.5)
    ax.axvline(x=0, color='#999999', linewidth=0.5)
    #ax.plot(V, [0.]*len(V), "-", color="#999999")
    ax.plot(V, [1e12*(F0 - gamma*v) for v in V], "-g", label=r"$F(0) - \gamma v$")
    ax.plot(V, [1e12*f for f in F], ".b", label=r"$F(v)$")

    ax.set_ylabel("force [pN]")
    ax.set_xlabel("velocity [m/s]")

    ax.legend(loc="best")
    #ax.grid()
    fig = plt.gcf()
    fig.set_size_inches((4,3))
    #nanopores.savefigs("howorka_v2f_2", FIGDIR)
    plt.show()

do_profile = False
if do_profile:
    # working 3D setup
    params = user_params(dim=3, Nmax=1.5e5, h=1., dnaqsdamp=0.25,
                         Qmol=-1., bulkcon=300.)
    # 2D setup
    #params = user_params(dim=2, Nmax=2e4, h=.5, dnaqsdamp=0.25,
    #                     Qmol=-1., bulkcon=300.)

    # along axis
    #Z = np.linspace(-6., 6., 42)
    #X = [[0.,0.,z] for z in Z]

    # at crosssection
    r0 = Howorka.params_geo3D.r0
    rMol = Howorka.default_geop.rMolecule
    eps = 1e-2
    R = r0 - rMol - eps
    Z = np.linspace(-R, R, 21)
    X = [[z,0.,0.] for z in Z]
    #X = [[0.,0.,0.]]
    print velocities(X, nproc=7, name="howorka_velo3D_2", **params)

do_plot = False
redo_plot = False

if do_plot:
    imax = user_params(imax=15)["imax"]
    if redo_plot:
        x = [0.2, 0., 0.]
        #x = [0., 0., 0.]
        setup = Howorka.Setup(x0=x, **params)
        f, v, dv = velocity_iteration(setup, imax)
        nanopores.save_stuff("velocity_iteration", f.tolist(), v.tolist(), dv.tolist())
    f, v, dv = nanopores.load_stuff("velocity_iteration")

    dim = params["dim"]
    plt.semilogy(range(1, imax+1), 1e12*np.sqrt(np.sum(np.array(f)**2, 1)),
                 "s-", label="net force on molecule")
    plt.ylabel("force [pN]")
    plt.xlabel("# iterations")
    plt.xlim(xmin=1, xmax=imax)
    plt.xticks(range(1,imax+1))
    plt.legend(loc="best")
    fig = plt.gcf()
    fig.set_size_inches((4,3))
    nanopores.savefigs("howorka_velocity_3D_z0", FIGDIR)