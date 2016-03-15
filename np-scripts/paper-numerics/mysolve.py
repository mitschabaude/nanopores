""" define solvers that record cumulative times needed for every loop """

from dolfin import *
from nanopores.tools.pdesystem import newtonsolve
from nanopores import *

def adaptive_pbpnps(geo, phys, cyl=False, frac=0.5, Nmax=1e4, Felref=None, Fdragref=None, Fpbref=None):
    LinearPB = LinearPBAxisymGoalOriented if cyl else LinearPBGoalOriented
    PNPStokes = PNPSAxisym if cyl else PNPS
    z = phys.dim - 1
    
    bV = phys.bV
    phys.bV = 0.
    goal = lambda v : phys.Fbare(v, z) #+ phys.CurrentPB(v)
    pb = LinearPB(geo, phys, goal=goal, ref=Fpbref)
    phys.bV = bV
    pb.maxcells = Nmax
    pb.marking_fraction = frac
    refined = True
    i = 0
    
    print "Number of cells:", pb.geo.mesh.num_cells()
        
    while refined:
        i += 1
        print "\n\nSolving PB."
        # solve pb
        pb.single_solve()
        pb.print_functionals()
        
        # define and solve pnps
        if i==1:
            pnps = PNPStokes(pb.geo, phys, v0=pb.solution)
        else:
            pnps.geo = pb.geo
            mesh = pb.geo.mesh
            for name, S in pnps.solvers.items():
                print "Adapting %s." % name
                S.adapt(mesh)
            functions = tuple(pnps.functions.values())
            for S in pnps.solvers.values():
                S.replace(functions,functions)
        print "\nSolving PNPS."
        dofs = sum(u.function_space().dim() for u in pnps.solutions())
        print "  Degrees of freedom: %d" % (dofs,)
        pnps.solve()
        #newton_iter = pnps.newton_solve()
        #print "  Newton iterations:", newton_iter
        print
        pnps.visualize("pore")
        fs = pnps.get_functionals()
        Fdrag = fs["Fp%d" %z] + fs["Fshear%d" %z]
        Fel = fs["Fbarevol%d" %z]
        F = Fdrag + Fel
        print "Fbare [pN]:", Fel
        print "Fdrag [pN]:", Fdrag
        print "F     [pN]:", F
        if Felref is not None:
            pb.save_estimate("Fel", abs((Fel-Felref)/Felref), N=dofs)
            pb.save_estimate("Fdrag", abs((Fdrag-Fdragref)/Fdragref), N=dofs)
            Fref = Felref + Fdragref
            pb.save_estimate("F", abs((F-Fref)/Fref), N=dofs)
        
        print "\nAdaptive refinement."
        (ind, err) = pb.estimate()
        pb.save_estimate("Fpb est", err, N=dofs)
        refined = pb.refine(ind)
        if not refined:
            print "Maximal number of cells reached."
        else:
            print "New total number of cells:", pb.geo.mesh.num_cells()

    return pb, pnps
    
def adaptive_pb(geo, phys, cyl=False, frac=0.5, Nmax=1e4, Fpbref=None, mesh2D=None, cheapest=False):
    LinearPB = LinearPBAxisymGoalOriented if cyl else LinearPBGoalOriented
    z = phys.dim - 1
    bV = phys.bV
    phys.bV = 0.
    goal = lambda v : phys.Fbare(v, z) #+ phys.CurrentPB(v)
    pb = LinearPB(geo, phys, goal=goal, ref=Fpbref)
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
        # solve pb
        pb.single_solve()
        pb.print_functionals()
        
        #plot(pb.geo.mesh)
        plot(pb.geo.submesh("membrane"))
        #plot(pb.geo.submesh("dna"))
        plot(pb.geo.submesh("pore"))
        #plot(pb.geo.submesh("molecule"))
        if mesh2D is not None:
            plot_cross(pb.functions["primal"], mesh2D, title="pb primal")
            plot_cross(pb.functions["dual"], mesh2D, title="pb dual")
        
        print "\nError estimation."
        (ind, err) = pb.estimate()
        print "\nMesh refinement."
        refined = pb.refine(ind)
        if not refined:
            print "Maximal number of cells reached."
        else:
            print "New total number of cells:", pb.geo.mesh.num_cells()
    return pb


def newton_solve(self, tol=None, damp=None, verbose=True):
    if tol is None: tol = self.tolnewton
    if damp is None: damp = self.newtondamp
    S = self.solvers.values()[0]    
    S.newtondamp = damp
    tcum = 0.
    Uold = self.solutions(deepcopy=True)
    
    for i in range(self.imax):
        tloop = Timer("loop")
        S.solve()
        if verbose:
            print '     Relative L2 Newton error:',S.relerror()
        if S.convergence(tol):
            if verbose:
                print "     Break loop because tolerance %s was reached." %tol
            converged = True
            break
        # cumulative time 
        tcum += tloop.stop()
        # calculate the error
        U = self.solutions(deepcopy=True)
        err = sum(errornorm(u, uold, "L2", degree_rise=0) for u, uold in zip(U, Uold)) / sum(norm(u, "L2") for u in U) 
        Uold = U
        self.save_estimate("err newton i", err, N=i+1)
        self.save_estimate("err newton time", err, N=tcum)
        continue
    else:
        if verbose: print "     Did not reach tolerance %s." %tol
        converged = False
    print "     Newton iterations:",i+1
        #print '     Relative L2 Newton error:',S.relerror()
    return i+1, converged
    

def hybrid_solve(self, tol=None, damp=None):
    if tol is None: tol = self.params["tolnewton"]
    if damp is None: damp = self.params["damp"]
    I = self.params["ipicard"]
    J = self.params["inewton"]
    nverbose = self.params["nverbose"]
    verbose = self.params["verbose"]
    times = {name : 0. for name in self.solvers}
    tcum = 0.
    Uold = self.solutions(deepcopy=True)

    for i in range(1, I+1):
        #v = self.functions["npm"]
        #plot(v)
        #self.visualize()
        if verbose:
            print "\n-- Fixed-Point Loop %d of max. %d" % (i, I)
        tloop = Timer("loop")
        for name, solver in self.solvers.items():
            if verbose:
                print "    Solving %s." % name
            t = Timer(name)
            if solver.is_linear:
                solver.solve()
                times[name] += t.stop()
            else:
                j, con = newtonsolve(solver, tol, damp, J, nverbose)
                times[name] += t.stop()
                if j==1 and con:
                    print "- Break at iteration %d because Newton stopped changing." %i
                    break
        else:
            # cumulative time 
            tcum += tloop.stop()
            # calculate the error
            U = self.solutions(deepcopy=True)
            err = sum(errornorm(u, uold, "L2", degree_rise=0) for u, uold in zip(U, Uold)) / sum(norm(u, "L2") for u in U) 
            Uold = U
            self.save_estimate("err hybrid i", err, N=i)
            self.save_estimate("err hybrid time", err, N=tcum)
            continue
        break
        
    Tt = sum(times.values())
    if verbose:
        print "\n CPU Time (solve): %.2f s" % Tt
        for name in self.solvers:
            print "  -) %s: %.2f s" %(name, times[name])
            
            
