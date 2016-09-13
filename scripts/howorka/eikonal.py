from dolfin import *
from forcefield2D import maybe_calculate
from nanopores.models import Howorka
import nanopores

def LinearPotential(rMol=0.5, relRange=1.25, U0=2.):
    d = relRange*rMol # range
    U0 = U0*nanopores.kT # height of barrier at r=rMol
    def U(r):
        if r > d:
            return 0.
        else:
            return U0*(r-d)/(rMol-d)
    return U

def boundary_force(mesh=None, rMolecule=0.5, **params):
    geo, phys = Howorka.setup2D(mesh=mesh, z0=None, **params)
    mesh = geo.submesh("fluid")
    geo, phys = Howorka.setup2D(mesh=mesh, z0=None, **params)
    
    mesh = geo.mesh
    V = FunctionSpace(mesh, "CG", 1)
    v = TestFunction(V)
    u = TrialFunction(V)
    f = Constant(1.0)
    y = Function(V)
    
    bc = geo.BC(V, Constant(0.), "dnab").bcs
    
    # Initialization problem to get good initial guess for nonlinear problem:
    F1 = inner(grad(u), grad(v))*dx - f*v*dx
    solve(lhs(F1)==rhs(F1), y, bc)
    
    # Stabilized Eikonal equation
    print "max cell size:", mesh.hmax()
    eps = Constant(mesh.hmax()/25)
    F = sqrt(inner(grad(y), grad(y)))*v*dx - f*v*dx + eps*inner(grad(y), grad(v))*dx
    # also works:
    #F = inner(grad(y), grad(y))*v*dx - f*v*dx + eps*inner(grad(y), grad(v))*dx
    solve(F==0, y, bc)
    
    print "Max distance:", y.vector().max()
    
    U = LinearPotential(rMol=rMolecule)
            
    class Potential(Expression):
        def eval(self, value, x):
            value[0] = U(y(x))
            
    U1 = Function(V)
    U1.interpolate(Potential())
    return -phys.grad(U1), U1

if __name__ == "__main__":
    F, Fel, Fdrag, mesh, params = maybe_calculate()
    F, U = boundary_force(mesh=mesh, **params)
    #F, U = boundary_force()
    plot(U, rescale=True, title="Potential")
    plot(F, rescale=True, title="Force")
    interactive()
