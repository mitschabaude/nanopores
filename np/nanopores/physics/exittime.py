from dolfin import *
from ..tools.pdesystem import GeneralLinearProblem
from ..tools.transientpde import *
from ..tools.illposed import IllposedLinearSolver, Functional

__all__ = ["ExitTimeProblem", "SurvivalProblem", "SuccessfulExit"]

class ExitTimeProblem(GeneralLinearProblem):
    k = 1
    
    method = dict(
        reuse = True,
        iterative = False,
        lusolver = ("superlu_dist" if has_lu_solver_method("superlu_dist") else "default"),
        luparams = dict(
            symmetric = False,
            same_nonzero_pattern = True,
            reuse_factorization = True,),
        ks = "bicgstab",
        kp = "amg",
        kparams = dict(
            maximum_iterations = 200,
            monitor_convergence = False,
            relative_tolerance = 1e-4,
            error_on_nonconvergence = False,
            preconditioner = dict(
                ilu = dict(fill_level = 1)))
    )

    @staticmethod
    def space(mesh):
        return FunctionSpace(mesh, 'CG', ExitTimeProblem.k)

    @staticmethod
    def forms(V, geo, phys, F):
        u = TrialFunction(V)
        v = TestFunction(V)
        dx = geo.dx("exittime")
        grad = phys.grad
        # TODO: for some reason, taking the pwconst causes conflict with F, results in u=NaN
        D = phys.DtargetBulk #geo.pwconst("Dtarget")
        mu = D/phys.kT
        
        J = -D*grad(v) + v*mu*F
        a = inner(J, grad(u))*dx # this is the dual of -div(J(u))
        L = Constant(-1.0)*v*dx
        
        return (a, L)
        
    @staticmethod
    def bcs(V, geo, exit={"exit"}):
        return [geo.BC(V, Constant(0.0), bou) for bou in exit]
        
        
class SurvivalProblem(GeneralLinearProblem):
    k = 1

    method = dict(
        reuse = True,
        iterative = False,
        lusolver = ("superlu_dist" if has_lu_solver_method("superlu_dist") else "default"),
        luparams = dict(
            symmetric = False,
            same_nonzero_pattern = True,
            reuse_factorization = True,),)
            
    @staticmethod
    def space(mesh):
        return FunctionSpace(mesh, 'CG', SurvivalProblem.k)

    @staticmethod
    def forms(V, geo, phys, u, F, dt):
        
        u1 = TrialFunction(V)
        v = TestFunction(V)
        dx = geo.dx("exittime")
        grad = phys.grad
        # TODO: for some reason, taking the pwconst causes conflict with F, results in u=NaN

        D = phys.DtargetBulk #geo.pwconst("Dtarget")
        mu = D/phys.kT
        J = -D*grad(v) + v*mu*F
        
        # (u1 - u)/dt + divJ*(th*u1 + (1-th)*u) = 0
        # u1 + dt*divJ*(u1*th) = u - dt*divJ*(u*(1-th))
        
        a = (u1*v - dt*inner(J, grad(u1)))*dx
        L = u*v*dx # + dt/2*inner(J, grad(u)))*dx
        
        return (a, L)
        
    @staticmethod
    def bcs(V, geo, exit={"exit"}):
        return [geo.BC(V, Constant(0.0), bou) for bou in exit]
        
    @staticmethod
    def initial_u(V):
        u = Function(V)
        # initial condition u(x,0) = 1
        u.vector()[:] = 1.
        return u
        
class SuccessfulExit(TransientLinearPDE):

    def __init__(self, geo=None, phys=None,
                 dt=None, badexit=set(), goodexit=set(), **problem_params):
        if dt is not None:
            self.dt = dt
            
        bothexit = badexit | goodexit
        
        badproblem = SurvivalProblem(geo=geo, phys=phys, dt=dt, exit=badexit, **problem_params)
        badsolver = IllposedLinearSolver(badproblem)
        ubad = badproblem.solution()
        
        bothproblem = SurvivalProblem(geo=geo, phys=phys, dt=dt, exit=bothexit, **problem_params)
        bothsolver = IllposedLinearSolver(bothproblem)
        uboth = bothproblem.solution()
        
        ugood = ubad - uboth

        self.geo = geo
        self.functions = {"bad": ubad, "both": uboth}
        self.solution = ugood
        self.solvers = {"bad": badsolver, "both": bothsolver}
        
        P = lambda z: ubad([0., 0., z]) - uboth([0., 0., z])
        zbtm = geo.params["zbtm"]
        ztop = geo.params["ztop"]
        
        self.plotter = TimeDependentPlotter(P, [zbtm, ztop, 200], dt)
        self.functionals = {}
        
    def visualize(self):
        self.plotter.plot(self.time[-1])
        
    def finish_plots(self):
        self.plotter.finish()
        self.plot_functionals("semilogx")
        
