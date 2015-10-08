""" Define Poisson problem """

from dolfin import *
from ..tools import AdaptableLinearProblem, LinearPDE, poisson_indicator
from .params_physical import *

parameters["allow_extrapolation"] = True
parameters["refinement_algorithm"] = "plaza_with_parent_facets"

__all__ = ["PoissonProblem", "Poisson", "PoissonProblemPureNeumannAxisym"]

class PoissonProblem(AdaptableLinearProblem):
    k = 1
    
    method = dict(
        reuse = True,
        iterative = True,
        lusolver = ("superlu_dist" if has_lu_solver_method("superlu_dist") else "default"),
        luparams = dict(
            symmetric = True,
            same_nonzero_pattern = True,
            reuse_factorization = True,),
        ks = "cg",
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
        return FunctionSpace(mesh, 'CG', PoissonProblem.k)

    @staticmethod
    def forms(V, geo, f):
        u = TrialFunction(V)
        v = TestFunction(V)
        dx = geo.dx()
        grad = geo.physics.grad
        lscale = geo.physics.lscale
        
        eps = geo.pwconst('permittivity')           
        a = inner(eps*grad(u), grad(v))*dx
        L = f*v*dx + lscale*geo.NeumannRHS(v, "surfcharge")
        
        return (a, L)
    
    def __init__(self, geo, phys=None, bcs=None, f=None, u=None):
        mesh = geo.mesh
        V = self.space(mesh)
        C0 = Constant(0.0)
        
        if f is None:
            f = C0
        if u is None:
            u = Function(V)
        if bcs is None:
            try:
                if phys.bV is None:
                    bcs = [geo.BC(V, C0, "ground")]
                else:
                    bcs = [geo.BC(V, phys.bV, "bV"),
                           geo.BC(V, C0, "ground")]
            except:
                warning("No boundary conditions have been assigned to %s" %type(self).__name__)

        a, L = self.forms(V, geo, f)

        AdaptableLinearProblem.__init__(self, a, L, u, bcs, geo.boundaries)

        
class PoissonProblemPureNeumannAxisym(PoissonProblem):
    
    method = PoissonProblem.method
    method["ks"] = "bicgstab"
    
    @staticmethod
    def space(mesh):
        V = FunctionSpace(mesh, "CG", PoissonProblemPureNeumannAxisym.k)
        R = FunctionSpace(mesh, "R", 0)
        return V * R
    
    @staticmethod
    def forms(V, geo, f):
        (u,b) = TrialFunctions(V)
        (v,c) = TestFunctions(V)
        dx = geo.dx()
        grad = geo.physics.grad
        lscale = geo.physics.lscale
        
        r = Expression("2*pi*x[0]")
        eps = geo.pwconst('permittivity')           
        a = inner(eps*grad(u), grad(v))*r*dx + (b*v + c*u)*r*dx
        
        L = f*v*r*dx + geo.linearRHS(v*r, "volcharge")
        L = L + lscale*geo.NeumannRHS(v*r, "surfcharge")
        
        return (a, L)
        
        
class Poisson(LinearPDE):
    def __init__(self, geo, **params):
        LinearPDE.__init__(self, geo, PoissonProblem, **params)
    def estimate(self):
        return poisson_indicator(self.geo, self.functions.values()[0])


