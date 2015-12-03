""" Define Poisson-Boltzmann problem """

from dolfin import *
from ..tools import AdaptableNonlinearProblem, LinearPDE, NonlinearPDE, poisson_indicator, GeneralNonlinearProblem, GeneralLinearProblem
from .params_physical import *

parameters["allow_extrapolation"] = True
parameters["refinement_algorithm"] = "plaza_with_parent_facets"

__all__ = ["NonstandardPB", "LinearNonstandardPB", "Poisson_"]

class NonstandardPBProblem(GeneralNonlinearProblem):
    k = 1
    
    method = dict(
        reuse = True,
        iterative = True,
        lusolver = ("superlu_dist" if has_lu_solver_method("superlu_dist") else "default"),
        luparams = dict(
            symmetric = True,
            same_nonzero_pattern = True,
            reuse_factorization = True,),
        ks = "bicgstab",
        kp = "amg",
        kparams = dict(
            maximum_iterations = 200,
            monitor_convergence = False,
            relative_tolerance = 1e-4,
            preconditioner = dict(
                ilu = dict(fill_level = 1)))
    )

    @staticmethod
    def space(mesh):
        return FunctionSpace(mesh, 'CG', NonstandardPBProblem.k)

    @staticmethod
    def forms(V, geo, phys, u):
        du = TrialFunction(V)
        v = TestFunction(V)
        dx = geo.dx()
         
        eps = geo.pwconst("permittivity")
        D = geo.pwconst("D")
        p0 = geo.pwconst("p0")
        n0 = geo.pwconst("n0")
        
        grad = phys.grad
        beta = 1./phys.UT
        
        #a = inner(eps*grad(U), grad(v))*dx + qq*beta*(n0 + p0)*U*v*dx
        #L = qq*D*v*dx + qq*(p0 - n0)*v*dx
        a = inner(eps*grad(du), grad(v))*dx + qq*beta*(n0*exp(beta*u) + p0*exp(-beta*u))*du*v*dx
        L = inner(eps*grad(u), grad(v))*dx - qq*(p0*exp(-beta*u) - n0*exp(beta*u))*v*dx - qq*D*v*dx
        return (a, L)
    
    @staticmethod
    def bcs(V, geo):
        return geo.pwconstBC(V, "bV")
        
    @staticmethod
    def bcs_homo(V, geo):
        return geo.pwconstBC(V, "bV", homogenize=True)
    
        
class LinearNonstandardPBProblem(GeneralLinearProblem):
    k = 1
    
    method = dict(
        reuse = True,
        iterative = True,
        lusolver = ("superlu_dist" if has_lu_solver_method("superlu_dist") else "default"),
        luparams = dict(
            symmetric = True,
            same_nonzero_pattern = True,
            reuse_factorization = True,),
        ks = "bicgstab",
        kp = "amg",
        kparams = dict(
            maximum_iterations = 200,
            monitor_convergence = False,
            relative_tolerance = 1e-4,
            preconditioner = dict(
                ilu = dict(fill_level = 1)))
    )

    @staticmethod
    def space(mesh):
        return FunctionSpace(mesh, 'CG', LinearNonstandardPBProblem.k)

    @staticmethod
    def forms(V, geo, phys, u):
        U = TrialFunction(V)
        v = TestFunction(V)
        dx = geo.dx()
         
        eps = geo.pwconst("permittivity")
        D = geo.pwconst("D")
        p0 = geo.pwconst("p0")
        n0 = geo.pwconst("n0")
        
        grad = phys.grad
        beta = 1./phys.UT
        
        a = inner(eps*grad(U), grad(v))*dx + qq*beta*(n0 + p0)*U*v*dx
        L = qq*D*v*dx + qq*(p0 - n0)*v*dx
        #a = inner(eps*grad(U), grad(v))*dx + qq*beta*(n0*exp(beta*u) + p0*exp(-beta*u))*U*v*dx
        #L = inner(eps*grad(u), grad(v))*dx + qq*(n0*exp(beta*u) - p0*exp(-beta*u))*v*dx - qq*D*v*dx
        return (a, L)
    
    @staticmethod
    def bcs(V, geo):
        return geo.pwconstBC(V, "bV")
        
class PoissonProblem(LinearNonstandardPBProblem):
    @staticmethod
    def forms(V, geo, phys, u):
        U = TrialFunction(V)
        v = TestFunction(V)
        dx = geo.dx()
         
        eps = geo.pwconst("permittivity")
        grad = phys.grad
        
        a = inner(eps*grad(U), grad(v))*dx
        L = Constant(0.)*v*dx
        return (a, L)
        

class NonstandardPB(NonlinearPDE):
    def __init__(self, geo, phys, **params):
        NonlinearPDE.__init__(self, geo, NonstandardPBProblem, phys=phys, **params)

class LinearNonstandardPB(LinearPDE):
    def __init__(self, geo, phys, **params):
        LinearPDE.__init__(self, geo, LinearNonstandardPBProblem, phys=phys, **params)
        
class Poisson_(LinearPDE):
    def __init__(self, geo, phys, **params):
        LinearPDE.__init__(self, geo, PoissonProblem, phys=phys, **params)

