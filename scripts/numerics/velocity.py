# (c) 2016 Gregor Mitscha-Baude
from nanopores.models import Howorka
from nanopores.tools import solvers

from nanopores import (LinearPBGoalOriented, LinearPBAxisymGoalOriented,
                       PNPSAxisym, PNPS)

def pbpnps(geo, phys, cyl=False, frac=0.5, Nmax=1e4, cheapest=False,
                      taylorhood=False):
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

def solve3D(geo, phys, **params):
    globals().update(params)
    if not stokesLU:
        StokesProblem.method["iterative"] = True
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
    return pbpnps(geo, phys, frac=frac3D, Nmax=Nmax3D, cheapest=cheapest)

# evaluate finite-size model for a number of z positions
def F_explicit3D(x, **params):
    import dolfin
    values = []
    for x0 in x:
        geo, phys = setup3D(x0=x0, **params)
        dolfin.plot(geo.submesh("solid"), key="b", title="solid mesh")
        pb, pnps = solve3D(geo, phys, **params)
        dolfin.plot(geo.submesh("solid"), key="b", title="solid mesh")
        values.append(pnps.forces())
    F, Fel, Fdrag = tuple(zip(*values))
    return F, Fel, Fdrag

print Howorka.F_explicit3D([[0.,0.,0.]], Nmax3D=4e5, taylorhood=False,
                           cheapest=True)

