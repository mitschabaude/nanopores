from dolfin import *
from forcefield2D import maybe_calculate
from nanopores.models import Howorka
import nanopores

F, Fel, Fdrag, mesh, params = maybe_calculate()
geo, phys = Howorka.setup2D(mesh=mesh, **params)
mesh = geo.submesh("fluid")
geo, phys = Howorka.setup2D(mesh=mesh, **params)

V = FunctionSpace(mesh, 'CG', 1)
v = TestFunction(V)
u = TrialFunction(V)
f = Constant(1.0)
y = Function(V)

bc = geo.BC(V, Constant(0.), "noslip").bcs

# Initialization problem to get good initial guess for nonlinear problem:
F1 = inner(grad(u), grad(v))*dx - f*v*dx
solve(lhs(F1)==rhs(F1), y, bc)

# Stabilized Eikonal equation
print "max cell size:", mesh.hmax()
eps = Constant(mesh.hmax()/25)
F = sqrt(inner(grad(y), grad(y)))*v*dx - f*v*dx + eps*inner(grad(y), grad(v))*dx
solve(F==0, y, bc)

print "Max distance:", y.vector().max()

# potential
rMol = 0.5
d = 1.25*rMol # range
U0 = 2.*nanopores.kT # height of barrier at r=rMol
def U(r):
    if r > d:
        return 0.
    else:
        return U0*(r-d)/(rMol-d)
        
def F(r):
    return 0. if r>d else U0/(rMol-d)
        
class Potential(Expression):
    def eval(self, value, x):
        value[0] = U(y(x))
        
U1 = Function(V)
U1.interpolate(Potential())

plot(-phys.grad(U1), rescale=True, title="Eikonal", interactive=True)
