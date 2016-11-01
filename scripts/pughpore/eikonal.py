from dolfin import *
from nanopores.tools import fields
import nanopores.geometries.pughpore as pughpore
import nanopores.models.pughpore as pughm
import nanopores as nano
fields.set_dir("/tmp/nanopores/")


#geop = nano.Params(pughpore.params)

def distance_boundary():
    up = nano.user_params(
        h=8.,
        R=30.,
        H=80.,
    )
    geop = nano.Params(pughpore.params)

    solverp = nano.Params(
        h = up.h,
        frac = 0.2,
        Nmax = 6e5,
        imax = 30,
        tol = 1e-2,
        cheapest = False,
    )
    physp = nano.Params()
    geo = pughpore.get_geo(solverp.h, **geop)
#    mesh = geo.submesh("fluid")
    phys = nano.Physics("default",geo,**physp)
#    geo, phys = Howorka.setup2D(mesh=mesh, z0=None, **params)
    
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
    eps = Constant(mesh.hmax()/25.)
    F = sqrt(inner(grad(y), grad(y)))*v*dx - f*v*dx + eps*inner(grad(y), grad(v))*dx
    # also works:
    #F = inner(grad(y), grad(y))*v*dx - f*v*dx + eps*inner(grad(y), grad(v))*dx
    solve(F==0, y, bc)
    
    print "Max distance:", y.vector().max()
    return y
    
    class Distance(Expression):
        def eval(self, value, x):
            value[0] = y(x)
            
    U = Function(V)
    U.interpolate(Distance())

    fields.save_functions("pugh_distance", {}, pugh_distance=U)
    fields.update()



    return U

if __name__ == "__main__":
    y = distance_boundary()
    setup=pughm.Setup()
    plotter=pughm.Plotter(setup)
    plotter.plot(y,title="Distance")
    interactive()
