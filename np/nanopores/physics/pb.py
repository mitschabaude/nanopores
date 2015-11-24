""" Define Poisson-Boltzmann problem """

from dolfin import *
from ..tools import AdaptableNonlinearProblem, NonlinearPDE, poisson_indicator, GeneralNonlinearProblem
from .params_physical import *

parameters["allow_extrapolation"] = True
parameters["refinement_algorithm"] = "plaza_with_parent_facets"

__all__ = ["PBProblem", "PB", "NonstandardPB"]

# TODO: cleanup

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
        ks = "cg",
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
        U = TrialFunction(V)
        v = TestFunction(V)
        dx = geo.dx()
         
        eps = geo.pwconst("permittivity")
        D = geo.pwconst("D")
        p0 = geo.pwconst("p0")
        n0 = geo.pwconst("n0")
        
        grad = phys.grad
        beta = 1./phys.UT
        
        a = inner(eps*grad(U), grad(v))*dx + qq*beta*(n0*exp(beta*u) + p0*exp(-beta*u))*U*v*dx
        L = inner(eps*grad(u), grad(v))*dx + qq*(n0*exp(beta*u) - p0*exp(-beta*u))*v*dx - qq*D*v*dx
        return (a, L)
    
    @staticmethod
    def bcs(V, geo):
        return geo.pwconstBC(V, "bV")
        
    @staticmethod
    def bcs_homo(V, geo):
        return geo.pwconstBC(V, "bV", homogenize=True)
 


class PBProblem(AdaptableNonlinearProblem):
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
            preconditioner = dict(
                ilu = dict(fill_level = 1)))
    )

    @staticmethod
    def space(mesh):
        return FunctionSpace(mesh, 'CG', PBProblem.k)

    @staticmethod
    def forms(V, geo, u0, phi):
        u = TrialFunction(V)
        v = TestFunction(V)
        dx = geo.dx()
        dx_ions = geo.dx("ions")
         
        eps = geo.pwconst('permittivity')
        phiF = geo.pwconst("fermi_level")
        c0 = geo.pwconst("ion_concentration")
        rho = geo.pwconst("permanent_charge")
        
        a = inner(eps*grad(u), grad(v))*dx + 2*cFarad*c0/UT*cosh((u0 - phiF)/UT)*u*v*dx
        L = inner(eps*grad(u0), grad(v))*dx + 2*cFarad*c0*sinh((u0 - phiF)/UT)*v*dx - rho*v*dx
        
        if geo.physicalboundary("charged"):
            L = L - phi*v('+')*geo.dS("charged")
        
        return (a, L)
    
    def __init__(self, geo, bcs=None, phi=None, u=None, **phys_params_in):
        mesh = geo.mesh
        V = self.space(mesh)
        phys_params.update(phys_params_in)
        C0 = Constant(0.0)
        
        if phi is None and geo.physicalboundary("charged"):
            area = assemble(Constant(1.0)*geo.dS("charged"))
            phi = Constant(phys_params["bcharge"]/area)
        if u is None:
            u = Function(V)
        if bcs is None:
            try:
                bcs = [geo.BC(V, phys_params["bV"], "bV"),
                       geo.BC(V, C0, "ground")]
            except:
                warning("No boundary conditions have been assigned to %s" %type(self).__name__)
                
        for bc in bcs:
            bc.apply(u.vector())
            #bc.homogenize() # <-- doesn't work!!!
        
        # TODO: decide on user interface    
        bcs = [geo.BC(V, C0, "bV"),
               geo.BC(V, C0, "ground")]
        
        a, L = self.forms(V, geo, u, phi)

        AdaptableNonlinearProblem.__init__(self, a, L, u, bcs, geo.boundaries)
        

class NonstandardPB(NonlinearPDE):
    def __init__(self, geo, phys, **params):
        NonlinearPDE.__init__(self, geo, NonstandardPBProblem, phys=phys, **params)
        
        
class PB(NonlinearPDE):
    def __init__(self, geo, **params):
        NonlinearPDE.__init__(self, geo, PBProblem, **params)
    def estimate(self):
        return poisson_indicator(self.geo, self.functions.values()[0])

