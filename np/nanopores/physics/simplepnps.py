""" Stripped-down, cleaner variants of PNPS allowing more general geometries """

from dolfin import *

from nanopores.tools import CoupledProblem, solvermethods, GeneralNonlinearProblem
from nanopores.physics.params_physical import *

#__all__ = ["SimplePNPS", "SimplePNPSAxisym"]
__all__ = ["SimplePNPProblem", "SimplePBProblem"]

# --- Problems ---

class SimplePNPProblem(GeneralNonlinearProblem):
    method = dict(solvermethods.bicgstab)
    method["iterative"] = False
    
    @staticmethod
    def space(mesh, k=1):
        V = FunctionSpace(mesh, 'CG', k)
        return MixedFunctionSpace((V, V, V))
        
    @staticmethod
    def initial_u(V, geo, phys, v0=None, vPB=None): # TODO incorporate initial guesses
        u = Function(V)
        u.interpolate(Constant((0.0, phys.bulkcon, phys.bulkcon)))
        return u

    @staticmethod
    def forms(V, geo, phys, u, ustokes=None, cyl=False):
        if not ustokes:
            dim = geo.mesh.topology().dim()
            ustokes = Constant(tuple(0. for i in range(dim)))
        
        dx = geo.dx()
        dx_ions = geo.dx("fluid")
        r2pi = Expression("2*pi*x[0]") if cyl else Constant(1.0)
        lscale = Constant(phys.lscale)
        grad = phys.grad

        eps = geo.pwconst("permittivity")
        Dp = geo.pwconst("Dp")
        Dm = geo.pwconst("Dm")
        kT = Constant(phys.kT)
        q = Constant(phys.qq)
        F = Constant(phys.cFarad)
        
        (v, cp, cm) = u.split()
        (w, dp, dm) = TestFunctions(V)
        
        apoisson = inner(eps*grad(v), grad(w))*r2pi*dx - F*(cp - cm)*w*r2pi*dx_ions
        aJm = inner(Dm*(grad(cm) - q/kT*cm*grad(v)) - cm*ustokes, grad(dm))*r2pi*dx_ions
        aJp = inner(Dp*(grad(cp) + q/kT*cp*grad(v)) - cp*ustokes, grad(dp))*r2pi*dx_ions
        
        Lqvol = geo.linearRHS(w*r2pi, "volcharge")
        Lqsurf = lscale*geo.NeumannRHS(w*r2pi, "surfcharge")
        
        L = apoisson + aJm + aJp - Lqvol - Lqsurf
        a = derivative(L, (v, cp, cm))

        return a, L
    
    @staticmethod
    def bcs(V, geo, phys):
        bcs = geo.pwconstBC(V.sub(0), "v0")
        bcs = bcs + geo.pwconstBC(V.sub(1), "c0")
        bcs = bcs + geo.pwconstBC(V.sub(2), "c0")
        return bcs
        
        
class SimplePBProblem(GeneralNonlinearProblem):
    method = dict(solvermethods.bicgstab)
    
    @staticmethod
    def space(mesh, k=1):
        return FunctionSpace(mesh, 'CG', k)

    @staticmethod
    def forms(V, geo, phys, u, cyl=False):
        dx = geo.dx()
        dx_ions = geo.dx("fluid")
        r2pi = Expression("2*pi*x[0]") if cyl else Constant(1.0)
        lscale = Constant(phys.lscale)
        grad = phys.grad

        eps = geo.pwconst("permittivity")
        UT = Constant(phys.UT)
        k = Constant(2*phys.cFarad*phys.bulkcon)
        
        w = TestFunction(V)
        
        apoisson = inner(eps*grad(u), grad(w))*r2pi*dx + k*sinh(u/UT)*w*r2pi*dx_ions
        Lqvol = geo.linearRHS(w*r2pi, "volcharge")
        Lqsurf = lscale*geo.NeumannRHS(w*r2pi, "surfcharge")
        
        L = apoisson - Lqvol - Lqsurf
        a = derivative(L, u)

        return a, L
    
    @staticmethod
    def bcs(V, geo, phys):
        # TODO: really only allow homogeneous BCs for PB?
        return geo.pwconstBC(V, "v0", homogenize=True)     

"""
class SimpleStokesProblemAxisym(StokesProblemAxisym):
    k = 1
    beta = 0.1

    @staticmethod
    def space(mesh):
        k = StokesProblemAxisymEqualOrder.k
        U = VectorFunctionSpace(mesh, 'CG', k)
        P = FunctionSpace(mesh, 'CG', k)
        return U*P

    @staticmethod
    def forms(W, geo, f):
        mesh = geo.mesh

        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)

        dx = geo.dx("fluid")
        r = Expression("x[0]")
        pi = 3.14159265359
        h = CellSize(mesh)
        delta = StokesProblemAxisymEqualOrder.beta*h**2

        # conservative formulation for correct BC, with added stabilization term
        a = (2*eta*inner(sym(grad(u)),sym(grad(v)))*r + 2*eta*u[0]*v[0]/r + \
            (div(v)*r+v[0])*p + q*(u[0] + div(u)*r))*2*pi*dx - \
            delta*inner(grad(p),grad(q))*r*2*pi*dx

        if f is None:
            f = Constant((0.0,0.0))

        L = 2*pi*inner(f,v - delta*grad(q))*r*dx
        return (a,L)
"""
