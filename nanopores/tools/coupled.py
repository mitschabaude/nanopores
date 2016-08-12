import dolfin
from collections import OrderedDict
from .utilities import _call
from .pdesystem import PDESystem, _pass, newtonsolve
from .illposed import IllposedLinearSolver, IllposedNonlinearSolver

__all__ = ["CoupledProblem", "CoupledSolver"]

class CoupledProblem(object):
    """ Automated creation of coupled problems out of single ones.
    Designed for problems that subclass General(Non)LinearProblem;
    will not immediately instantiate Problems but first use their static methods.

    FIXME: One should also be able to make CoupledProblems out of CoupledProblems,
        to create nested fixed point loops automatically.
        to accomplish this, the CoupledProblem would have automatically provide a
        .space, which can be a accessed from the class *before* instantiation
        -- completely unreasonable in this case.
        -> We should split the Problem initialisation API into two stages:
           -) one where only information are gathered
           -) a second where the Function is created
        Maybe even the whole Problem concept is bad... the core classes could
        be Equation, Function (for which an Equation can be solved) and Solver
    
    Expected usage:
        problems = OrderedDict([
            ("pnp", PNPProblem),
            ("stokes", StokesProblem)])
        
        def couple_pnp(ustokes):
            return dict(ustokes = ustokes.sub(0))
        
        def couple_stokes(upnp, phys):
            v, cp, cm = upnp.split()
            f = -phys.cFarad*(cp - cm)*grad(v)
            return dict(f = f)
        
        couplers = dict(
            pnp = couple_pnp,
            stokes = couple_stokes
        )
        
        coupled = CoupledProblem(problems, couplers, geo, phys, **params)
        solver = CoupledSolver(coupled)
        for n in range(5):
            solver.solve()
    
    The CoupledProblem will use the given problem names to create internal names for the solutions, e.g.
    upnp, ustokes in this case, which will be fed through the couplers to provide arguments for the
    the problem's form definitions. Everything contained in the **params may be used in the couplers.
    E.g. in the above situation the StokesProblem forms must only expect a variable "f" (and thus not know
    anything about its coupling to PNP):
    
    class StokesProblem(GeneralLinearProblem)
        @staticmethod
        def forms(..., f):
            L = inner(f, v)
            ...
            
    and CoupledProblem will make sure that f actually becomes -cFarad*(cp - cm)*grad(v).
    The use of OrderedDict in the example makes sure that PNPProblem will be solved first in every
    iteration. If this is unimportant, problems may be simply a dict.
    """

    def __init__(self, problems, couplers, geo, phys=None, **params):
        # problems = dictionary of Problem classes
        # keys are important! (see class docstring)
        
        self.solutions = OrderedDict()
        self.oldsolutions = OrderedDict()
        self.problems = OrderedDict()
        self.geo = geo
        
        params.update(geo=geo, phys=phys)
        self.params = params
        mesh = geo.mesh
        
        for name, Problem in problems.items():
            # get space
            V = _call(Problem.space, dict(params, mesh=mesh))
            
            # create solution function
            if hasattr(Problem, "initial_u"):
                u = _call(Problem.initial_u, dict(params, V=V))
            else:
                u = dolfin.Function(V)
            # and backup function for solution in last iteration
            uold = dolfin.Function(V)
                
            self.solutions[name] = u
            self.oldsolutions[name] = uold
            params.update({"u" + name: u})
            
        # now we have all solutions in hand and can actually create the problems
        for name, Problem in problems.items():
            u = self.solutions[name]
            uold = self.oldsolutions[name]
            
            # obtain parameters coupled to other problems
            problem_specific_params = dict(params, u=u, uold=uold)
            problem_specific_params.update(
                _call(couplers[name], problem_specific_params))
            
            # instantiate the problem
            self.problems[name] = Problem(**problem_specific_params)
            
    def update_uold(self):
        # useful to e.g. change timestep and reassemble matrices
        # this assumes that coupled parameters are *not* changed
        for name in self.problems:
            uold = self.oldsolutions[name]
            u = self.solutions[name]
            #uold.vector()[:] = u.vector()[:] # <-- does not work after adapting
            uold.assign(u.copy(deepcopy=True))
            
        
    def update_forms(self, **new_params):
        # useful to e.g. change timestep and reassemble matrices
        # this assumes that coupled parameters are *not* changed
        for name, problem in self.problems.items():
            problem.update_forms(**new_params)
        
class CoupledSolver(PDESystem):
    params = dict(inewton = 10, ipicard = 10,
        tolnewton = 1e-4, damp = 1., verbose=True, nverbose=False)

    def __init__(self, coupled, goals=[], **solverparams):
                
        self.solvers = OrderedDict()
        for name, problem in coupled.problems.items():
            if problem.is_linear:
                solver = IllposedLinearSolver(problem)
                solver.is_linear = True
            else:
                solver = IllposedNonlinearSolver(problem)
                solver.is_linear = False
            self.solvers[name] = solver
            
        self.geo = coupled.geo
        self.functions = coupled.solutions
        self.coupled = coupled
        self.problems = coupled.problems
        self.functionals = {}
        self.add_functionals(goals)
        self.params.update(solverparams)
    
    # TODO: clean up -- single_solve() should rely on fixedpoint()
    # TODO: explore possibility to choose newton tol adaptively    
    def single_solve(self, tol=None, damp=None, inside_loop=_pass):
        if tol is None: tol = self.params["tolnewton"]
        if damp is None: damp = self.params["damp"]
        I = self.params["ipicard"]
        J = self.params["inewton"]
        nverbose = self.params["nverbose"]
        verbose = self.params["verbose"]
        times = {name : 0. for name in self.solvers}
    
        for i in range(1, I+1):
            if verbose:
                print "\n-- Fixed-Point Loop %d of max. %d" % (i, I)
            for name, solver in self.solvers.items():
                if verbose:
                    print "    Solving %s." % name
                t = dolfin.Timer(name)
                if solver.is_linear:
                    solver.solve()
                    times[name] += t.stop()
                else:
                    j, con = newtonsolve(solver, tol, damp, J, nverbose, lambda: inside_loop(self))
                    times[name] += t.stop()
                    if j==1 and con:
                        print "- Break at iteration %d because Newton stopped changing." %i
                        break
            else: 
                inside_loop(self)
                self.coupled.update_uold()
                continue
            break
        Tt = sum(times.values())
        if verbose:
            print "\n CPU Time (solve): %.2f s" % Tt
            for name in self.solvers:
                print "  -) %s: %.2f s" %(name, times[name])
                   
    def fixedpoint(self, tol=None, damp=None):
        if tol is None: tol = self.params["tolnewton"]
        if damp is None: damp = self.params["damp"]
        imax = self.params["ipicard"]
        inewton = self.params["inewton"]
        nverbose = self.params["nverbose"]
        verbose = self.params["verbose"]     
        ntol = tol              
                       
        times = {name : 0. for name in self.solvers}
        tcum = 0.
        U = self.coupled.solutions
        Uold = self.coupled.oldsolutions
        self.converged = False
    
        for i in range(1, imax+1):
            if verbose:
                print "\n-- Fixed-Point Loop %d of max. %d" % (i, imax)
            tloop = dolfin.Timer("loop")
            for name, solver in self.solvers.items():
                if verbose:
                    print "    Solving %s." % name
                t = dolfin.Timer(name)
                if solver.is_linear:
                    solver.solve()
                else:
                    for _ in newtoniteration(
                        solver, ntol, damp, inewton, nverbose):
                        pass
                times[name] += t.stop()        
            # cumulative time 
            tcum += tloop.stop()
            
            yield i
            
            # calculate the error
            errors = [(name, error(U[name], Uold[name])) for name in U]
            if verbose:
                for item in errors: print "    error %s: %s" % item
            err = sum(err for _, err in errors)/len(errors)
            
            # check for stopping
            if err < tol:
                self.converged = True
                if verbose:
                    print "- Break at iteration %d because err %.3g < tol %.3g" %(i, err, tol)
                break
            
            self.save_estimate("err hybrid i", err, N=i)
            self.save_estimate("err hybrid time", err, N=tcum)
            self.coupled.update_uold()
            
        Tt = sum(times.values())
        if verbose:
            print "\n CPU Time (solve): %.2f s" % Tt
            for name in self.solvers:
                print "  -) %s: %.2f s" %(name, times[name])
                
    def generic_fixedpoint(self, tol=None, damp=None):
        # does timing, errors, stopping, but relies on user to implement solve
        if tol is None: tol = self.params["tolnewton"]
        if damp is None: damp = self.params["damp"]
        imax = self.params["ipicard"]
        verbose = self.params["verbose"]              
                       
        tcum = 0.
        U = self.coupled.solutions
        Uold = self.coupled.oldsolutions
        self.converged = False
    
        for i in range(1, imax+1):
            if verbose:
                print "\n-- Fixed-Point Loop %d of max. %d" % (i, imax)
            tloop = dolfin.Timer("loop")
            # solve here
            yield i
            # cumulative time 
            tcum += tloop.stop()
            
            # calculate the error
            errors = [(name, error(U[name], Uold[name])) for name in U]
            if verbose:
                for item in errors: print "    error %s: %s" % item
            err = sum(err for _, err in errors)/len(errors)
            
            # check for stopping
            if err < tol:
                self.converged = True
                if verbose:
                    print "- Break at iteration %d because err %.3g < tol %.3g" %(i, err, tol)
                break
            
            self.save_estimate("err hybrid i", err, N=i)
            self.save_estimate("err hybrid time", err, N=tcum)
            self.coupled.update_uold()

                
# for fixed point error criterion
def error(u, uold):
    norm = dolfin.norm(u, "L2")
    if norm==0.:
        return 1.
    return dolfin.errornorm(u, uold, "L2", degree_rise=0)/norm    
                
TOL = 1e-4
IMAX = 10

def newtoniteration(solver, tol=TOL, damp=1., imax=IMAX, verbose=True):
    "minimal newton iteration on one newton solver"
    solver.newtondamp = damp
    solver.converged = False
    
    for i in range(imax):
        solver.solve()
        if solver.convergence(tol):
            solver.converged = True
        yield solver
        if solver.converged:
            break              
                            
def plainfixedpoint(solvers, tol=TOL, ntol=TOL, imax=IMAX, inewton=1, 
                        damp=1., verbose=True, nverbose=False):                         
    "minimal fixed point iteration on coupled solver."
    for i in range(1, imax+1):
        for solver in solvers.values():
            if solver.is_linear:
                solver.solve()
            else:
                for _ in newtoniteration(solver, ntol, damp, inewton, nverbose):
                    pass
        yield solvers
