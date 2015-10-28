from dolfin import *
from ..tools.pdesystem import GeneralLinearProblem

__all__ = ["ExitTimeProblem"]

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
        D = phys.Dtarget # TODO: pwconst
        mu = D/phys.kT
        
        a = -inner(D*grad(u), grad(v))*dx + inner(mu*F, grad(u))*v*dx
        L = Constant(-1.0)*v*dx
        
        return (a, L)
        
    @staticmethod
    def bcs(V, geo):
        return [geo.BC(V, Constant(0.0), "exit")]

