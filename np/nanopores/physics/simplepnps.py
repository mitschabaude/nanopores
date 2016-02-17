""" Stripped-down, cleaner variants of PNPS allowing more general geometries """

from dolfin import *

from nanopores.tools import CoupledProblem, solvermethods, GeneralNonlinearProblem
from nanopores.physics.params_physical import *

#__all__ = ["SimplePNPS", "SimplePNPSAxisym"]
__all__ = ["SimplePNPProblem"]

# --- Problems ---

class SimplePNPProblem(GeneralNonlinearProblem):
    method = solvermethods.bicgstab
    method["iterative"] = False
    
    @staticmethod
    def space(mesh, k=1):
        V = FunctionSpace(mesh, 'CG', k)
        return MixedFunctionSpace((V, V, V))
        
    @staticmethod
    def initial_u(V, geo, phys, v0=None, vPB=None): # TODO
        u = Function(V)
        u.interpolate(Constant((0.0, phys.bulkcon, phys.bulkcon)))
        return u

    @staticmethod
    def forms(V, geo, phys, u, ustokes=None, axisymmetric=False):
        # TODO: automatic Jacobian
        (v, cp, cm) = TrialFunctions(V)
        (vv, dp, dm) = TestFunctions(V)
        
        (vold, cpold, cmold) = u.split()

        if not ustokes:
            dim = geo.mesh.topology().dim()
            ustokes = Constant(tuple(0. for i in range(dim)))
        uold = ustokes
        
        dx = geo.dx()
        dx_ions = geo.dx("fluid")
        r2pi = Expression("2*pi*x[0]") if axisymmetric else Constant(1.0)

        eps = geo.pwconst("permittivity")
        Dp = geo.pwconst("Dp")
        Dm = geo.pwconst("Dm")
        kT = phys.kT
        q = phys.qq
        F = phys.cFarad
        
        apoisson = inner(eps*grad(v),grad(vv))*r2pi*dx - F*(cp - cm)*vv*r2pi*dx_ions
        aJm = inner(Dm*(grad(cm) - q/kT*(cm*grad(vold) + cmold*grad(v))) - cm*uold, grad(dm))*r2pi*dx_ions
        aJp = inner(Dp*(grad(cp) + q/kT*(cp*grad(vold) + cpold*grad(v))) - cp*uold, grad(dp))*r2pi*dx_ions
        a = apoisson + aJm + aJp

        Lpoisson = inner(eps*grad(vold),grad(vv))*r2pi*dx - F*(cpold - cmold)*vv*r2pi*dx_ions
        LJm = inner(Dm*(grad(cmold) - q/kT*cmold*grad(vold)) - cmold*uold, grad(dm))*r2pi*dx_ions
        LJp = inner(Dp*(grad(cpold) + q/kT*cpold*grad(vold)) - cpold*uold, grad(dp))*r2pi*dx_ions
        Lqvol = geo.linearRHS(vv*r2pi, "volcharge")
        Lqsurf = geo.NeumannRHS(vv*r2pi, "surfcharge")
        Lq = Lqvol + Lqsurf

        L = Lpoisson + LJm + LJp - Lq
        
        return a, L
    
    @staticmethod
    def bcs(V, geo, phys):
        
        c0 = Constant(phys.bulkcon)
        bcs = geo.pwconstBC(V.sub(0), "v0")
        bcs = bcs + geo.pwconstBC(V.sub(1), "c0")
        bcs = bcs + geo.pwconstBC(V.sub(2), "c0")
        return bcs

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
