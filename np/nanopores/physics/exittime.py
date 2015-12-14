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
        kp = "ilu",
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
            reuse_factorization = True,),
        ks = "bicgstab",
        kp = "ilu",
        kparams = dict(
            maximum_iterations = 1000,
            monitor_convergence = False,
            relative_tolerance = 1e-4,
            error_on_nonconvergence = False,
            preconditioner = dict(
                ilu = dict(fill_level = 1)))
    )
            
    @staticmethod
    def space(mesh):
        return FunctionSpace(mesh, 'CG', SurvivalProblem.k)

    @staticmethod
    def forms(V, geo, phys, u, F, dt=None, steady=False):
        
        u1 = TrialFunction(V)
        v = TestFunction(V)
        dx = geo.dx("exittime")
        grad = phys.grad
        # TODO: for some reason, taking the pwconst causes conflict with F, results in u=NaN

        D = phys.DtargetBulk #geo.pwconst("Dtarget")
        mu = D/phys.kT
        J = -D*grad(v) + v*mu*F
        
        if steady or dt is None:
            a = inner(J, grad(u1))*dx
            L = Constant(0.)*v*dx
        else: # transient case
            # backward euler: (u1 - u)/dt + divJ*(u1) = 0
            a = (u1*v - dt*inner(J, grad(u1)))*dx
            L = u*v*dx
        
        return (a, L)
        
    @staticmethod
    def bcs(V, geo, goodexit=set(), badexit={"exit"}):
        return ([geo.BC(V, Constant(0.0), bou) for bou in badexit] +
                [geo.BC(V, Constant(1.0), bou) for bou in goodexit])
        
    @staticmethod
    def initial_u(V, u0=0.):
        u = Function(V)
        # initial condition u(x,0) = u0
        # u0 = 1 for "survival" type problem, u0 = 0 for "death"
        u.vector()[:] = u0
        return u

class SuccessfulExit(TransientLinearPDE):

    def __init__(self, geo=None, phys=None, dt=None, **problem_params):
        TransientLinearPDE.__init__(self, SurvivalProblem, geo=geo,
            phys=phys, dt=dt, **problem_params)
        
        p = self.solution
        pz = lambda z: p([0., 0., z])
        zbtm = geo.params["zporebtm"]
        ztop = geo.params["ztop"]
        
        self.plotter = TimeDependentPlotter(pz, [zbtm, ztop, 200], dt)
        
    def visualize(self):
        self.plotter.plot(self.time[-1])
        
    def finish_plots(self):
        self.plotter.finish()
        self.plot_functionals("semilogx")
        
