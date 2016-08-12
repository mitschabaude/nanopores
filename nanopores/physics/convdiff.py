"standard convection-diffusion/transport equation with a given force field"
"u_t = div(-D*grad(u) + D/kT*F*u) + f"

from dolfin import *
from nanopores.tools.pdesystem import GeneralLinearProblem, LinearPDE
from nanopores.tools.transientpde import *
from nanopores.tools import solvermethods

__all__ = ["ConvectionDiffusionProblem", "ConvectionDiffusion"]
                
class ConvectionDiffusionProblem(GeneralLinearProblem):
    method = dict(solvermethods.direct_reuse)
            
    @staticmethod
    def space(mesh, k=1, steady=False):
        V = FunctionSpace(mesh, 'CG', k)
        if steady:
            R = FunctionSpace(mesh, 'Real', 0)
            return V*R
        else:
            return V
        
    @staticmethod
    def initial_u(V, u0=None, steady=False):
        u = Function(V)
        if steady and u0 is not None:            
            W = V.sub(0).collapse()
            w = interpolate(u0, W)
            R = V.sub(1).collapse()
            assign(u, [w, Function(R)])
        elif u0 is not None:
            u.interpolate(u0)
        return u

    @staticmethod
    def forms(V, geo, phys, u, F, f=None, dt=None, steady=False, cyl=False):
        
        u1 = TrialFunction(V)
        v = TestFunction(V)
        if steady:
            u1, c = split(u1)
            v, d = split(v)
            u0, _ = split(u)
            
        r2pi = Expression("2*pi*x[0]") if cyl else Constant(1.0)
        dx = geo.dx("fluid")
        grad = phys.grad
        #lscale = Constant(phys.lscale)
        #n = FacetNormal(geo.mesh)

        D = geo.pwconst("Dtarget")
        #D = Constant(phys.DtargetBulk)
        kT = Constant(phys.kT)
        # J = -D*grad(u) + D/kT*F*u
        def J(u):
            return -D*grad(u) + D/kT*F*u
        
        if f is None: # injection function
            f = Constant(0.)
        
        if steady or dt is None:
            # -div(J) = f mit Neumann constraint
            a = inner(J(u1), grad(v))*r2pi*dx + (c*v + u1*d)*r2pi*dx
            L = f*v*r2pi*dx + u0*d*r2pi*dx
        else:
            dt = Constant(dt)
            # u_t = -div(J(u)) - f
            # backward euler:
            # (u1 - u)/dt = -div(J(u1)) - f
            # u1 + dt*div(J(u1)) = u - dt*f
            a = (u1*v - dt*inner(J(u1), grad(v)))*r2pi*dx
            L = (u*v - dt*f*v)*r2pi*dx
            
            #aNoBC = dt*lscale*inner(J(u1), n*v)*r2pi*geo.ds("lowerb")
        
        return (a, L)
        
    @staticmethod
    def bcs(V, geo, bc={}, steady=False):
        if steady:
            return geo.pwBC(V.sub(0), "c0", value=bc)
        else:
            return geo.pwBC(V, "c0", value=bc)
        
class ConvectionDiffusion(TransientLinearPDE):

    def __init__(self, geo=None, phys=None, dt=None, **problem_params):
        TransientLinearPDE.__init__(self, ConvectionDiffusionProblem, geo=geo,
            phys=phys, dt=dt, **problem_params)
            
class ConvectionDiffusionSteady(LinearPDE):

    def __init__(self, geo=None, phys=None, **problem_params):
        LinearPDE.__init__(self, geo, ConvectionDiffusionProblem, 
            phys=phys, steady=True, **problem_params)
        
