""" General-purpose PDE System class """

import numpy
from dolfin import *
from .illposed import *
from .errorest import *
from .utilities import _call
from collections import OrderedDict

parameters["refinement_algorithm"] = "plaza_with_parent_facets"

__all__ = ["PDESystem", "LinearPDE", "NonlinearPDE", "GoalAdaptivePDE",
           "GeneralLinearProblem", "GeneralNonlinearProblem",
           "solve_pde", "solve_problem", "PDEfromProblem"]
                 
_pass = lambda *args, **kwargs : None

class PDESystem(object):
    imax = 100
    maxcells = 10000
    marking_fraction = 0.8
    uniform_refinement = False

    def __init__(self, geo=None, solvers={}, functions={}, functionals={}):
        self.geo = geo
        self.functions = functions
        self.solvers = solvers
        self.functionals = functionals

    def solve(self, refinement=False, verbose=True, inside_loop=_pass):

        if verbose:
            print("Number of cells:",self.geo.mesh.num_cells())
        if self.geo.mesh.num_cells() > self.maxcells:
            refinement = False

        for i in range(self.imax):
            if verbose:
                print('\n- Loop ' +str(i+1) + ' of max.', self.imax)
                print("  Degrees of freedom: %d" % (sum(u.function_space().dim() for u in self.solutions()),))

            self.single_solve() #inside_loop=inside_loop)
            if verbose:
                self.print_functionals()
            if inside_loop is not None:
                inside_loop(self)
            #plot(self.solvers.values()[0].problem.solution(), interactive=True)
            #self.visualize()
            if refinement:
                (ind,err) = self.estimate()
                self.save_estimate("est", err)
                if verbose:
                    print("Error estimate (H1):",err)
                refined = self.refine(ind)
                if not refined:
                    if verbose:
                        print('Maximal number of cells reached',  \
                               ' \n  ==> no more refinement \n')
                    break
                elif verbose:
                    print("New total number of cells:",self.geo.mesh.num_cells())
            else:
                break

    def save_estimate(self, string, err, N=None):
        if not hasattr(self, "estimators"):
            self.estimators = {}
        if string not in self.estimators:
            self.estimators[string] = Estimator(string)
        if N is None:
            N = self.geo.mesh.num_cells()
        self.estimators[string] += N, err

    def estimate_uniform(self):
        "the trivial indicator"
        # for uniform refinement 
        # TODO: does not work with doerfler marking
        return None, 1.
    
    def estimate_zz(self):
        """ simple zz indicator, estimator """
        # just to have a basic estimation tool when nothing is defined
        u = self.solutions()[0]
        mesh = self.geo.mesh
        ind, err = zz_indicator(u)
        return ind, err
        
    estimate = estimate_zz #uniform

    def single_solve(self, **other):
        for S in list(self.solvers.values()): S.solve()

    def refine(self, ind):
        mesh0 = self.geo.mesh
        mesh = self.refine_mesh(ind)

        #plot(mesh0)
        #plot(mesh, interactive=True)
        #self.adapt(mesh)

        if mesh.num_cells() > self.maxcells:
            #self.geo.mesh = mesh0
            return False

        self.adapt(mesh)
        return True

    def refine_mesh(self, ind):
        mesh = self.geo.mesh
        
        # MARK
        markers = CellFunction("bool", mesh, True)
        if not self.uniform_refinement and not self.marking_fraction == 1.:
            #tic()
            indicators = CellFunction("double", mesh)
            # ind is a DG0 Function
            dofmap = ind.function_space().dofmap()
            cell_to_dof = numpy.array([dofmap.cell_dofs(i)[0] for i in range(mesh.num_cells())])
            indicators.array()[:] = ind.vector()[cell_to_dof]
            
            # This version is about 30x slower than the one above:
            # (and only about 0.5x faster than the refine() step
            #  which should dominate this function)
            #for c in cells(mesh):
            #    indicators[c] = ind(c.midpoint())
            #print "TIME DG0 -> CellFunction: %s s" % (toc(),)
            #tic()
            dorfler_mark(markers, indicators, self.marking_fraction)
            #print "TIME Marking: %s s" % (toc(),)
        #tic()

        # REFINE
        mesh = refine(mesh, markers)
        #print "TIME Refinement: %s s" % (toc(),)
        #mesh = adapt(mesh, markers)
        return mesh

    def adapt(self, mesh):
        self.geo.adapt(mesh)
        
        for name, S in list(self.solvers.items()):
            #print "Adapting %s." % name
            S.adapt(mesh)

        #for S in self.solvers.values():
        #    S.adapt(mesh)

        functions = tuple(self.functions.values())
        for S in list(self.solvers.values()):
            S.replace(functions,functions)

        for J in list(self.functionals.values()):
            if isinstance(J,list):
                for j in J:
                    j.adapt(mesh)
                    j.replace(functions,functions)
            else:
                J.adapt(mesh)
                J.replace(functions,functions)

    def rebuild(self, mesh):
        """ Assumes geometry to have geo.rebuild """
        functionals = self.functionals

        self.geo.rebuild(mesh)
        self.__init__(self.geo)

        for Jstr,J in list(self.functionals.items()):
            J.values = functionals[Jstr].values

    def visualize(self, subdomain=None):
        sol = {x: self.solutions(x, deepcopy=True) for x in self.functions}
        mesh = self.geo.mesh
        on = ""

        if subdomain:
            on = " on " + subdomain
            mesh = self.geo.submesh(subdomain)

        plot(mesh, title="final mesh"+on)
        for x in sol:
            for i, f in enumerate(sol[x]):
                if subdomain:
                    adaptfunction(f, mesh, assign=True)
                plot(f, title = ("%s-%i" % (x, i)) +on)
        interactive()
        
    def add_functionals(self, functionals):
        # the input functionals here is a list of functions of the solutions (flat tuple) and geo
        # these should each return a dict(name of functional = ufl form)
        U = self.solutions()
        geo = self.geo
        for f in functionals:
            self.functionals.update({key: Functional(F) for key, F in list(f(U, geo).items())})

    def print_functionals(self):
        Jdir = self.functionals
        for Jstr in sorted(self.functionals):
            J = Jdir[Jstr]
            if isinstance(J,list):
                for ii in range(len(J)):
                    print(("%s[%i]: " %(Jstr,ii)) + str(J[ii].evaluate()))
            else:
                print(("%s: " %Jstr) + str(J.evaluate()))

    def get_functional(self, Jstr):
        return self.functionals[Jstr].evaluate()
    
    def get_functionals(self, Jlist=None):
        if Jlist is None:
            Jlist = sorted(self.functionals.keys())
        Jdic = self.functionals
        return OrderedDict([(J, Jdic[J].evaluate()) for J in Jlist])
        """
        functionals = {}
        for Jstr in Jlist:
            J = Jdir[Jstr]
            if isinstance(J,list):
                for ii in range(len(J)):
                    functionals["%s[%i]" %(Jstr,ii)] = J[ii].evaluate()
            else:
                functionals[Jstr] = J.evaluate()
        return functionals
        """
        
    def evaluate(self, f):
        return f(self.solutions())

    def print_results(self, names=None):
        if not names:
            self.print_functionals()

    def solutions(self, string=None, deepcopy=False):
        if string:
            f = self.functions[string]
            if f.function_space().num_sub_spaces() > 0:
                return f.split(deepcopy=deepcopy)
            else:
                return (f,)
        t = ()
        for x in self.functions:
            t = t + self.solutions(x, deepcopy=deepcopy)
        return t

    def dofs(self):
        return sum(u.function_space().dim() for u in self.solutions())

    def save_mesh(self, mesh_name=None):
        geo_name = self.geo.parameter("name")
        from nanopores import DATADIR
        DIR = "%s/%s/mesh" %(DATADIR, geo_name)
        if not mesh_name:
            mesh_name = "last_adapted_mesh"

        meshfile = File("%s/%s.xml" %(DIR, mesh_name))
        meshfile << self.geo.mesh
        N = str(self.geo.mesh.num_cells())
        meshfile = File("%s/adapted/mesh_%s.xml" %(DIR, N))
        meshfile << self.geo.mesh
        return DIR
        
def newtonsolve(S, tol=1e-4, damp=1., imax=10, verbose=True, inside_loop=_pass): 
    S.newtondamp = damp
    for i in range(imax):
        S.solve()
        #plot(self.solution) # for debugging
        inside_loop()
        if verbose:
            print('     Relative L2 Newton error:',S.relerror())
        if S.convergence(tol):
            if verbose:
                print("     Break loop because tolerance %s was reached." %tol)
            converged = True
            break
    else:
        if verbose: print("     Did not reach tolerance %s." %tol)
        converged = False
    print("     Newton iterations:",i+1)
        #print '     Relative L2 Newton error:',S.relerror()
    return i+1, converged
    
    
class LinearPDE(PDESystem):
    ''' simple interface for single linear PDE '''
    def __init__(self, geo, ProblemClass, *problem_args, **problem_params):
        problem = ProblemClass(geo, *problem_args, **problem_params)
        solver = IllposedLinearSolver(problem)

        self.geo = geo
        self.functions = {ProblemClass.__name__: problem.solution()}
        self.solution = problem.solution()
        self.problem = problem
        self.solvers = {ProblemClass.__name__: solver}
        self.functionals = {}

class NonlinearPDE(PDESystem):
    ''' simple interface for single nonlinear PDE and Newton method '''
    tolnewton = 1e-4
    newtondamp = 1.

    def __init__(self, geo, ProblemClass, **problem_params):
        problem = ProblemClass(geo, **problem_params)
        solver = IllposedNonlinearSolver(problem)

        self.geo = geo
        self.functions = {ProblemClass.__name__: problem.solution()}
        self.solution = problem.solution()
        self.problem = problem
        self.solvers = {ProblemClass.__name__: solver}
        self.functionals = {}
        
    def single_solve(self, tol=None, damp=None, imax=None, verbose=True, inside_loop=_pass):
        if not tol: tol = self.tolnewton
        if not damp: damp = self.newtondamp
        if not imax: imax = self.imax
        S = list(self.solvers.values())[0]
        return newtonsolve(S, tol, damp, imax, verbose, lambda: inside_loop(self))
        
      
def solve_pde(Problem, geo=None, phys=None, refinement=False, imax = 20, maxcells=1e4,
        marking_fraction=0.8, tolnewton=1e-2, newtondamp=1., iterative=None, visualize=False,
        inside_loop=_pass, goals=(), **params):
    """ very simple interface for quick tests """
    solverparams = dict(imax=imax, maxcells=maxcells, marking_fraction=marking_fraction,
        tolnewton=tolnewton, newtondamp=newtondamp)
    if iterative is not None:
        Problem.method["iterative"] = iterative # TODO shouldn't be necessary to change class attributes
        
    PDEClass = LinearPDE if Problem.is_linear else NonlinearPDE
    pde = PDEClass(geo, Problem, phys=phys, **params)
    for key in solverparams:
        setattr(pde, key, solverparams[key])
    pde.add_functionals(goals)
        
    t = Timer("solve")
    pde.solve(refinement=refinement, inside_loop=inside_loop)
    #pde.single_solve(inside_loop=inside_loop)
    print("CPU time (solve): %s [s]" % (t.stop(),))
    
    if visualize:
        pde.visualize()
    return pde
    

class PDEfromProblem(LinearPDE, NonlinearPDE):

    def __init__(self, problem, geo):
        import types
        if problem.is_linear:
            solver = IllposedLinearSolver(problem)
            self.single_solve = types.MethodType(LinearPDE.single_solve, self)
        else:
            solver = IllposedNonlinearSolver(problem)
            self.single_solve = types.MethodType(NonlinearPDE.single_solve, self)
        self.geo = geo
        self.functions = {type(problem).__name__: problem.solution()}
        self.solution = problem.solution()
        self.problem = problem
        self.solvers = {type(problem).__name__: solver}
        self.functionals = {}

        
def solve_problem(problem, geo, imax = 20, maxcells=1e4,
        marking_fraction=0.8, tolnewton=1e-2, newtondamp=1., iterative=None, visualize=False,
        goals=(), **solve_params):
    "simple interface for quick tests. like solve_pde, but takes instantiated problem; useful for customized problems"
    solverparams = dict(imax=imax, maxcells=maxcells, marking_fraction=marking_fraction,
        tolnewton=tolnewton, newtondamp=newtondamp)
    if iterative is not None:
        problem.method["iterative"] = iterative
        
    pde = PDEfromProblem(problem, geo)
    for key in solverparams:
        setattr(pde, key, solverparams[key])
    pde.add_functionals(goals)
        
    t = Timer("solve")
    pde.solve(**solve_params)
    print("CPU time (solve): %s [s]" % (t.stop(),))
    
    if visualize:
        pde.visualize()
    return pde
    

class GoalAdaptivePDE(PDESystem):
    ''' simple interface for PDE solver with goal-oriented adaptivity '''
    def __init__(self, geo, phys, Problem, goal):
        # create two problems: the primal one with forms a, L and the *dual* one with
        # a_dual(u,v) := a(v,u)
        # L_dual(v) := goal(v)
        primal = Problem(geo, phys)
        dual = Problem(geo, phys)

        # now dual is just another instance of the primal problem, but we modify it:
        aT = adjoint(dual.a)
        dual.a = aT
        # now let v be the TestFunction of the adjoint problem
        v = aT.arguments()[0]
        # create new RHS, the goal
        L = goal(v)
        dual.L = L

        solver = IllposedLinearSolver(primal)
        dualsolver = IllposedLinearSolver(dual)

        # in the end, we are mainly interested in the goal functional
        # evaluated at the solution of the primal problem
        # so this should be a Functional as well
        u = primal.solution()
        goal_f = Functional(goal(u))

        self.geo = geo
        self.phys = phys
        self.functions = {"primal": u, "dual":dual.solution()}
        self.solution = u
        self.solvers = {"primal": solver, "dual":dualsolver}
        self.functionals = {"goal": goal_f}
        
    
class GeneralLinearProblem(AdaptableLinearProblem):
    is_linear = True

    def __init__(self, geo, phys=None, u=None, bcs=None, **params):
        
        mesh = geo.mesh
        V = _call(self.space, dict(params, mesh=mesh))
        params.update(geo=geo, phys=phys, V=V)
        
        if not u:
            if hasattr(self, "initial_u"):
                u = _call(self.initial_u, params)
            else:
                u = Function(V)
            
        params.update(u=u)
        self.params = params
        
        if not bcs:
            bcs = _call(self.bcs, params)
        
        a, L = _call(self.forms, params)
        AdaptableLinearProblem.__init__(self, a, L, u, bcs, geo.boundaries)
        
    def update_forms(self, **new_params):
        # useful to e.g. change timestep and reassemble matrix
    
        self.params.update(new_params)
        a, L = _call(self.forms, self.params)
        self.a = a
        self.L = L

        
class GeneralNonlinearProblem(AdaptableNonlinearProblem):
    is_linear = False

    def __init__(self, geo, phys=None, u=None, bcs=None, **params):
        
        mesh = geo.mesh
        V = _call(self.space, dict(params, mesh=mesh))
        params.update(geo=geo, phys=phys, V=V)
        
        if not u:
            if hasattr(self, "initial_u"):
                u = _call(self.initial_u, params)
            else:
                u = Function(V)
            
        params.update(u=u)
        self.params = params
        
        if not bcs:
            bcs = _call(self.bcs, params)
        
        a, L = _call(self.forms, params)
        AdaptableNonlinearProblem.__init__(self, a, L, u, bcs, geo.boundaries)

