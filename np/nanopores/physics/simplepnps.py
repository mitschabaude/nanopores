""" Stripped-down, cleaner variants of PNPS allowing more general geometries """

from dolfin import *

from nanopores.tools import CoupledProblem, solvermethods, GeneralNonlinearProblem, GeneralLinearProblem
from nanopores.physics.params_physical import *

#__all__ = ["SimplePNPS", "SimplePNPSAxisym"]
#__all__ = ["SimplePNPProblem", "SimplePBProblem"]

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
        if ustokes is None:
            dim = geo.mesh.topology().dim()
            ustokes = Constant(tuple(0. for i in range(dim)))
        
        dx = geo.dx()
        dx_ions = geo.dx("fluid")
        n = FacetNormal(geo.mesh)
        r2pi = Expression("2*pi*x[0]", degree=1) if cyl else Constant(1.0)
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
        
        Jm = -Dm*(grad(cm) - q/kT*cm*grad(v)) + cm*ustokes
        Jp = -Dp*(grad(cp) + q/kT*cp*grad(v)) + cp*ustokes
        
        apoisson = inner(eps*grad(v), grad(w))*r2pi*dx - F*(cp - cm)*w*r2pi*dx_ions
        aJm = inner(Jm, grad(dm))*r2pi*dx_ions
        aJp = inner(Jp, grad(dp))*r2pi*dx_ions
        
        # TODO: investigate "no bcs" further. in the test problem, they don't work as expected
        aNoBCp = -jump(lscale*Jp*dp*r2pi, n)*geo.dS("nocbc") - lscale*inner(Jp, n*dp)*r2pi*geo.ds("nocbc")
        aNoBCm = -jump(lscale*Jm*dm*r2pi, n)*geo.dS("nocbc") - lscale*inner(Jm, n*dm)*r2pi*geo.ds("nocbc")
        
        Lqvol = geo.linearRHS(w*r2pi, "volcharge")
        Lqsurf = lscale*geo.NeumannRHS(w*r2pi, "surfcharge")
        LJm = lscale*geo.NeumannRHS(dm*r2pi, "cmflux")
        LJp = lscale*geo.NeumannRHS(dp*r2pi, "cpflux")
        
        L = apoisson + aJm + aJp + aNoBCp + aNoBCm - Lqvol - Lqsurf - LJm - LJp
        a = derivative(L, (v, cp, cm))

        return a, L
    
    @staticmethod
    def bcs(V, geo, phys):
        bcs = geo.pwBC(V.sub(0), "v0")
        bcs += geo.pwBC(V.sub(1), "cp0")
        bcs += geo.pwBC(V.sub(2), "cm0")
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
        
        
class SimplePoissonProblem(GeneralLinearProblem):
    method = dict(solvermethods.poisson)
    
    @staticmethod
    def space(mesh, k=1):
        return FunctionSpace(mesh, 'CG', k)

    @staticmethod
    def forms(V, geo, phys, u, cyl=False):
        dx = geo.dx()
        r2pi = Expression("2*pi*x[0]") if cyl else Constant(1.0)
        lscale = Constant(phys.lscale)
        grad = phys.grad
        eps = geo.pwconst("permittivity")
        
        w = TestFunction(V)
        apoisson = inner(eps*grad(u), grad(w))*r2pi*dx
        Lqvol = geo.linearRHS(w*r2pi, "volcharge")
        Lqsurf = lscale*geo.NeumannRHS(w*r2pi, "surfcharge")
        L = apoisson - Lqvol - Lqsurf
        a = derivative(L, u)
        return a, L
    
    @staticmethod
    def bcs(V, geo, phys):
        return geo.pwconstBC(V, "v0")  


class SimpleStokesProblem(GeneralLinearProblem):
    "stabilized equal-order formulation; consistent for k=1"
    method = dict(solvermethods.stokes)

    @staticmethod
    def space(mesh, ku=1, kp=1):
        U = VectorFunctionSpace(mesh, 'CG', ku)
        P = FunctionSpace(mesh, 'CG', kp)
        return U*P

    @staticmethod
    def forms(V, geo, phys, f=None, cyl=False, beta=0.01, conservative=True):
        # beta = stabilization parameter, TODO: better lower in 2D?
        mesh = geo.mesh
        if f is None:
            dim = geo.mesh.topology().dim()
            f = Constant(tuple(0. for i in range(dim)))

        (u, p) = TrialFunctions(V)
        (v, q) = TestFunctions(V)

        grad = phys.grad
        div = phys.div
        lscale = phys.lscale
        
        dx = geo.dx("fluid")
        r = Expression("x[0]")
        pi2 = Constant(2.*pi)
        h = CellSize(mesh)
        delta = Constant(beta/lscale**2)*h**2
        eta = Constant(phys.eta)
        
        def eps(u): return Constant(2.)*sym(grad(u))

        # conservative formulation for correct BC, with added stabilization term
        if cyl:
            a = (eta*inner(eps(u), eps(v))*r + Constant(2.)*eta*u[0]*v[0]/r + \
                (div(v)*r+v[0])*p + q*(u[0] + div(u)*r))*pi2*dx - \
                delta*inner(grad(p), grad(q))*r*pi2*dx
            L = inner(f, v - delta*grad(q))*r*pi2*dx
        else:
            a = (eta*inner(eps(u), eps(v)) + div(v)*p + q*div(u))*dx \
                 - delta*inner(grad(p), grad(q))*dx
            L = inner(f, v - delta*grad(q))*dx
            
        if not conservative and cyl:
            a = (eta*inner(grad(u), grad(v))*r + eta*u[0]*v[0]/r - \
                inner(v, grad(p))*r + q*(u[0] + div(u)*r))*pi2*dx - \
                delta*inner(grad(p), grad(q))*r*pi2*dx
            L = inner(f, v - delta*grad(q))*r*pi2*dx
        if not conservative and not cyl:
            a = (eta*inner(grad(u), grad(v)) - inner(v, grad(p)) + q*div(u))*dx \
                - delta*inner(grad(p), grad(q))*dx
            L = inner(f, v - delta*grad(q))*dx
        
            #FIXME
            #p = 2*inner(sym(grad(u)), sym(grad(v)))*dx + lscale*inner(p, q)*dx
        # TODO: include preconditioning form in some way
        return a, L
        
    @staticmethod
    def bcs(V, geo):
        return geo.pwBC(V.sub(0), "noslip") + geo.pwBC(V.sub(1), "pressure")
        
        


