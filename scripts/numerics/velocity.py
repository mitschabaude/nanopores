# (c) 2016 Gregor Mitscha-Baude
import dolfin, os
from nanopores.models import Howorka
from nanopores.tools import solvers
import nanopores.tools.fields as fields
from nanopores import (LinearPBGoalOriented, LinearPBAxisymGoalOriented,
                       PNPSAxisym, PNPS, StokesProblem, PNPProblem, HOME)
DATADIR = os.path.join(HOME, "Dropbox", "nanopores", "fields")
fields.set_dir(DATADIR)

def pbpnps(geo, phys, cyl=False, frac=0.5, Nmax=1e4, cheapest=False,
                      taylorhood=False, stokesLU=False, **kwargs):
    if not stokesLU and not cyl:
        StokesProblem.method["iterative"] = True
    else:
        StokesProblem.method["iterative"] = False
    PNPProblem.method["iterative"] = False
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

# evaluate finite-size model for a number of z positions
@solvers.cache_forcefield("howorka_nonzero_u")
def F(x, dim=3, UMol=None, **params):
    values = []
    setup = Howorka.setup3D if dim==3 else Howorka.setup2D
    cyl = dim==2
    for x0 in x:
        geo, phys = setup(x0=x0, **params)
        if UMol is not None:
            phys.update(UMol=UMol)
        #dolfin.plot(geo.submesh("solid"), key="b", title="solid mesh")
        pb, pnps = pbpnps(geo, phys, cyl=cyl, **params)
        #dolfin.plot(geo.submesh("solid"), key="b", title="solid mesh")
        values.append(pnps.forces())
        pnps.visualize("pore")
    F, Fel, Fdrag = tuple(zip(*values))
    return dict(F=F, Fel=Fel, Fdrag=Fdrag)

print F([[0.,0.,0.]], Nmax=2e4, UMol=(0.,0.1,-0.05), dim=3, dnaqsdamp=0.5,
        taylorhood=False, cheapest=True, cache=False, h3D=2.,
        stokesLU=True)

