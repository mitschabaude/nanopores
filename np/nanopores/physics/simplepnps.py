""" Stripped-down, cleaner variants of PNPS allowing more general geometries """
import numpy
from dolfin import *
from collections import OrderedDict
from nanopores.tools import (CoupledProblem, solvermethods,
    GeneralNonlinearProblem, GeneralLinearProblem, CoupledSolver)

__all__ = ["SimplePNPProblem", "SimplePBProblem", "SimpleStokesProblem",
           "SimplePoissonProblem",
           "PNPSHybrid", "PNPSFixedPoint", "PNPFixedPoint",
           "PNPSFixedPointbV", "PNPFixedPointNonlinear"]

# --- Problems ---

class SimplePNPProblem(GeneralNonlinearProblem):
    method = dict(solvermethods.bicgstab)
    method["iterative"] = False
    
    @staticmethod
    def space(mesh, k=1):
        V = FunctionSpace(mesh, 'CG', k)
        return MixedFunctionSpace((V, V, V))
        
    @staticmethod
    def initial_u(V, geo, phys, v0=None): # TODO incorporate initial guesses
        u = Function(V)
        if v0 is None:
            u.interpolate(Constant((0.0, phys.bulkcon, phys.bulkcon)))
        else:
            W = V.sub(0).collapse()
            v = interpolate(v0, W)
            c0 = phys.bulkcon
            cp = Function(W)
            cm = Function(W)
            cp.vector()[:] = c0*numpy.exp(-v.vector()[:]/phys.UT)
            cm.vector()[:] = c0*numpy.exp(v.vector()[:]/phys.UT)
            assign(u, [v, cp, cm])
        return u

    @staticmethod
    def forms(V, geo, phys, u, ustokes=None, cyl=False):
        if ustokes is None:
            dim = geo.mesh.topology().dim()
            ustokes = Constant(tuple(0. for i in range(dim)))
        
        dx = geo.dx()
        dx_ions = geo.dx("fluid")
        n = FacetNormal(geo.mesh)
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
        return geo.pwBC(V.sub(0), "v0") + \
               geo.pwBC(V.sub(1), "cp0") + \
               geo.pwBC(V.sub(2), "cm0")
               
        
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
    def forms(V, geo, phys, f=None, dxf=None, cyl=False):
        dx = geo.dx()
        r2pi = Expression("2*pi*x[0]") if cyl else Constant(1.0)
        lscale = Constant(phys.lscale)
        grad = phys.grad
        eps = geo.pwconst("permittivity")
        
        v = TrialFunction(V)
        w = TestFunction(V)
        
        a = inner(eps*grad(v), grad(w))*r2pi*dx
        Lqvol = geo.linearRHS(w*r2pi, "volcharge")
        Lqsurf = lscale*geo.NeumannRHS(w*r2pi, "surfcharge")
        L = Lqvol + Lqsurf
        if f is not None:
            L = L + f*w*r2pi*dxf
        return a, L
    
    @staticmethod
    def bcs(V, geo, phys):
        return geo.pwBC(V, "v0")
        
class LinearSGPoissonProblem(GeneralLinearProblem):
    "Linearized Scharfetter-Gummel-type Poisson problem for fixed point PNP"
    method = dict(solvermethods.bicgstab)
    
    @staticmethod
    def space(mesh, k=1):
        return FunctionSpace(mesh, 'CG', k)

    @staticmethod
    def forms(V, geo, phys, u, cp, cm, dx_ions, cyl=False):
        dx = geo.dx()
        r2pi = Expression("2*pi*x[0]") if cyl else Constant(1.0)
        lscale = Constant(phys.lscale)
        grad = phys.grad
        eps = geo.pwconst("permittivity")
        UT = Constant(phys.UT)
        F = Constant(phys.cFarad)
        
        v = TrialFunction(V)
        w = TestFunction(V)
        
        a = inner(eps*grad(v), grad(w))*r2pi*dx + F*v/UT*(cp + cm)*w*r2pi*dx_ions
        
        Lqions = F*((cp - cm) + u/UT*(cp + cm))*w*r2pi*dx_ions
        Lqvol = geo.linearRHS(w*r2pi, "volcharge")
        Lqsurf = lscale*geo.NeumannRHS(w*r2pi, "surfcharge")
        L = Lqions + Lqvol + Lqsurf
        return a, L
    
    @staticmethod
    def bcs(V, geo, phys):
        return geo.pwBC(V, "v0")

class SGPoissonProblem(GeneralNonlinearProblem):
    "Scharfetter-Gummel-type Poisson problem for fixed point PNP"
    method = dict(solvermethods.bicgstab)
    
    @staticmethod
    def space(mesh, k=1):
        return FunctionSpace(mesh, 'CG', k)

    @staticmethod
    def forms(V, geo, phys, u, uold, cp, cm, dx_ions, cyl=False):
        dx = geo.dx()
        r2pi = Expression("2*pi*x[0]") if cyl else Constant(1.0)
        lscale = Constant(phys.lscale)
        grad = phys.grad
        eps = geo.pwconst("permittivity")
        UT = Constant(phys.UT)
        F = Constant(phys.cFarad)

        w = TestFunction(V)
        
        apoisson = inner(eps*grad(u), grad(w))*r2pi*dx \
            - F*(exp(-(u-uold)/UT)*cp - exp((u-uold)/UT)*cm)*w*r2pi*dx_ions

        Lqvol = geo.linearRHS(w*r2pi, "volcharge")
        Lqsurf = lscale*geo.NeumannRHS(w*r2pi, "surfcharge")
        
        L = apoisson - (Lqvol + Lqsurf)
        a = derivative(L, u)
        return a, L
    
    @staticmethod
    def bcs(V, geo, phys):
        return geo.pwBC(V, "v0")        
        
class SimpleNernstPlanckProblem(GeneralLinearProblem):
    method = dict(solvermethods.bicgstab)
    
    @staticmethod
    def space(mesh, k=1):
        return FunctionSpace(mesh, "CG", k)
        
    @staticmethod
    def initial_u(V, phys):
        u = Function(V)
        u.interpolate(Constant(phys.bulkcon))
        return u

    @staticmethod
    def forms(V, geo, phys, z, E, D=None, ustokes=None, cyl=False):
        if ustokes is None:
            dim = phys.dim
            ustokes = Constant(tuple(0. for i in range(dim)))
        if D is None:
            D = geo.pwconst("D")
            
        dx = geo.dx("fluid")
        r2pi = Expression("2*pi*x[0]") if cyl else Constant(1.0)
        grad = phys.grad
        kT = Constant(phys.kT)
        q = Constant(phys.qq)
        
        c = TrialFunction(V)
        d = TestFunction(V)
        
        J = -D*grad(c) + z*q*D/kT*E*c + c*ustokes
        a = inner(J, grad(d))*r2pi*dx
        L = Constant(0.)*d*dx
        return a, L
            
    @staticmethod
    def bcs(V, geo, phys):
        return geo.pwBC(V, "c0")  


class SimpleStokesProblem(GeneralLinearProblem):
    "stabilized equal-order formulation; consistent for k=1"
    method = dict(solvermethods.stokes)

    @staticmethod
    def space(mesh, ku=1, kp=1):
        U = VectorFunctionSpace(mesh, 'CG', ku)
        P = FunctionSpace(mesh, 'CG', kp)
        return U*P

    @staticmethod
    def forms(V, geo, phys, f=None, cyl=False, beta=.01, conservative=True):
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
        r = Expression("x[0]/L", L=Constant(lscale))
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
            
        # optional non-conservative formulation with neumann BC n*grad(u) = 0
        if not conservative and cyl:
            a = (eta*inner(grad(u), grad(v))*r + eta*u[0]*v[0]/r - \
                inner(v, grad(p))*r + q*(u[0] + div(u)*r))*pi2*dx - \
                delta*inner(grad(p), grad(q))*r*pi2*dx
            L = inner(f, v - delta*grad(q))*r*pi2*dx
        if not conservative and not cyl:
            a = (eta*inner(grad(u), grad(v)) - inner(v, grad(p)) + q*div(u))*dx \
                - delta*inner(grad(p), grad(q))*dx
            L = inner(f, v - delta*grad(q))*dx
        
        # TODO: be able to include preconditioning form
        # p = 2*inner(sym(grad(u)), sym(grad(v)))*dx + lscale*inner(p, q)*dx
        return a, L
        
    @staticmethod
    def bcs(V, geo):
        return geo.pwBC(V.sub(0), "noslip") + geo.pwBC(V.sub(1), "pressure")
        
# --- hybrid solver ---

class PNPSHybrid(CoupledSolver):
    
    def __init__(self, geo, phys, goals=[], iterative=False, **params):
        problems = OrderedDict([
            ("pnp", SimplePNPProblem),
            ("stokes", SimpleStokesProblem),])

        def couple_pnp(ustokes):
            return dict(ustokes = ustokes.sub(0))
        def couple_stokes(upnp, phys):
            v, cp, cm = upnp.split()
            f = -phys.cFarad*(cp - cm)*phys.grad(v)
            return dict(f = f)
        couplers = dict(pnp=couple_pnp, stokes=couple_stokes)
    
        problem = CoupledProblem(problems, couplers, geo, phys, **params)
        problem.problems["pnp"].method["iterative"] = iterative
        problem.problems["stokes"].method["iterative"] = iterative
        CoupledSolver.__init__(self, problem, goals, **params)
    
# --- PNP fixed-point solvers ---
    
class PNPFixedPoint(CoupledSolver):

    def __init__(self, geo, phys, goals=[], iterative=False, **params):
        problems = OrderedDict([
            ("poisson", LinearSGPoissonProblem),
            ("npp", SimpleNernstPlanckProblem),
            ("npm", SimpleNernstPlanckProblem),
            ])

        def couple_poisson(unpp, unpm, geo):
            return dict(cp=unpp, cm=unpm, dx_ions=geo.dx("fluid"))   
        def couple_npp(upoisson, geo, phys):
            E = -phys.grad(upoisson)
            D = geo.pwconst("Dp")
            return dict(z=1., E=E, D=D)
        def couple_npm(upoisson, geo, phys):
            E = -phys.grad(upoisson)
            D = geo.pwconst("Dm")
            return dict(z=-1., E=E, D=D)
            
        couplers = dict(poisson=couple_poisson, npp=couple_npp, npm=couple_npm)
    
        problem = CoupledProblem(problems, couplers, geo, phys, **params)
        for name in problems:
            problem.problems[name].method["iterative"] = iterative
        CoupledSolver.__init__(self, problem, goals, **params)
        
class PNPFixedPointNonlinear(CoupledSolver):

    def __init__(self, geo, phys, goals=[], iterative=False, **params):
        problems = OrderedDict([
            ("poisson", SGPoissonProblem),
            ("npp", SimpleNernstPlanckProblem),
            ("npm", SimpleNernstPlanckProblem),
            ])

        def couple_poisson(upoisson, unpp, unpm, geo):
            return dict(v0=upoisson, cp=unpp, cm=unpm,
                        dx_ions=geo.dx("fluid"))   
        def couple_npp(upoisson, geo, phys):
            E = -phys.grad(upoisson)
            D = geo.pwconst("Dp")
            return dict(z=1., E=E, D=D)
        def couple_npm(upoisson, geo, phys):
            E = -phys.grad(upoisson)
            D = geo.pwconst("Dm")
            return dict(z=-1., E=E, D=D)
            
        couplers = dict(poisson=couple_poisson, npp=couple_npp, npm=couple_npm)
    
        problem = CoupledProblem(problems, couplers, geo, phys, **params)
        for name in problems:
            problem.problems[name].method["iterative"] = iterative
        CoupledSolver.__init__(self, problem, goals, **params)

class PNPFixedPointNaive(CoupledSolver):

    def __init__(self, geo, phys, goals=[], iterative=False, **params):
        problems = OrderedDict([
            ("poisson", SimplePoissonProblem),
            ("npp", SimpleNernstPlanckProblem),
            ("npm", SimpleNernstPlanckProblem),
            ])

        def couple_poisson(unpp, unpm, geo, phys):
            f = phys.cFarad*(unpp - unpm)
            dxf = geo.dx("fluid")
            return dict(f=f, dxf=dxf)   
        def couple_npp(upoisson, geo, phys):
            E = -phys.grad(upoisson)
            D = geo.pwconst("Dp")
            return dict(z=1., E=E, D=D)
        def couple_npm(upoisson, geo, phys):
            E = -phys.grad(upoisson)
            D = geo.pwconst("Dm")
            return dict(z=-1., E=E, D=D)
            
        couplers = dict(poisson=couple_poisson, npp=couple_npp, npm=couple_npm)
    
        problem = CoupledProblem(problems, couplers, geo, phys, **params)
        for name in problems:
            problem.problems[name].method["iterative"] = iterative
        CoupledSolver.__init__(self, problem, goals, **params)
        
# --- PNPS fixed-point solvers ---        

class PNPSFixedPoint(CoupledSolver):

    def __init__(self, geo, phys, goals=[], iterative=False, **params):
        problems = OrderedDict([
            ("poisson", LinearSGPoissonProblem),
            ("npp", SimpleNernstPlanckProblem),
            ("npm", SimpleNernstPlanckProblem),
            ("stokes", SimpleStokesProblem),
        ])

        def couple_poisson(unpp, unpm, geo):
            return dict(cp=unpp, cm=unpm, dx_ions=geo.dx("fluid"))     
        def couple_npp(upoisson, ustokes, geo, phys):
            E = -phys.grad(upoisson)
            D = geo.pwconst("Dp")
            return dict(z=1., E=E, D=D, ustokes=ustokes.sub(0))
        def couple_npm(upoisson, ustokes, geo, phys):
            E = -phys.grad(upoisson)
            D = geo.pwconst("Dm")
            return dict(z=-1., E=E, D=D, ustokes=ustokes.sub(0))
        def couple_stokes(upoisson, unpp, unpm, phys):
            v, cp, cm = upoisson, unpp, unpm
            f = -phys.cFarad*(cp - cm)*phys.grad(v)
            return dict(f = f)
            
        couplers = dict(poisson=couple_poisson, npp=couple_npp,
                        npm=couple_npm, stokes=couple_stokes)
    
        problem = CoupledProblem(problems, couplers, geo, phys, **params)
        for name in problems:
            problem.problems[name].method["iterative"] = iterative
        CoupledSolver.__init__(self, problem, goals, **params)
        
    def solve_pnp(self):
        for pde in "poisson", "npp", "npm":
            self.solvers[pde].solve()
            
    def solve_stokes(self):
        self.solvers["stokes"].solve()
        
    def fixedpoint(self, ipnp=2):
        for i in self.generic_fixedpoint():
            self.solve_pnp()
            if i > ipnp:
                self.solve_stokes()
            yield i
          
# bVscheme needs access to voltage bias = phys.bV
# this would hamper generality of PNPSFixedPoint class
import math

class PNPSFixedPointbV(PNPSFixedPoint):
    "voltage is slowly increased and stokes solved only afterwards"
    
    def fixedpoint(self, bVstep=0.025):
        bV = self.coupled.params["phys"].bV
        
        idamp = math.ceil(abs(bV)/bVstep)
        damping = 1./idamp if idamp != 0 else 1.
        ipnp = idamp + 1

        for i in self.generic_fixedpoint():
            if i <= idamp:
                self.solvers["poisson"].damp_bcs(damping*min(i, idamp))
                
            self.solve_pnp()
            if i > ipnp:
                self.solve_stokes()
            yield i

