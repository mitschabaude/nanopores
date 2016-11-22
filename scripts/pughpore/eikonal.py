"solve eikonal equation to get distance to boundary"
"TODO: implement also in 2D"

from dolfin import *
from nanopores.tools import fields
import nanopores.geometries.pughpore as pughpore
import nanopores.models.pughpore as pughm
import nanopores as nano
import folders

up = nano.user_params(h=8.0)

def distance_boundary():
    h = up.h
    geo = pughpore.get_geo(h)
    y = distance_boundary_from_geo(geo)
    print "Max distance:", y.vector().max()

    if not fields.exists("pugh_distance", h=h):
        fields.save_functions("pugh_distance", dict(h=h), y=y)
        fields.update()
    return y

def distance_boundary_from_geo(geo):
    mesh = geo.mesh
    V = FunctionSpace(mesh, "CG", 1)
    v = TestFunction(V)
    u = TrialFunction(V)
    f = Constant(1.0)
    y = Function(V)

    bc = geo.BC(V, Constant(0.), "dnab").bcs

    # Initialization problem to get good initial guess for nonlinear problem:
    F1 = inner(grad(u), grad(v))*dx - f*v*dx
    problem = LinearVariationalProblem(lhs(F1), rhs(F1), y, bc)
    solver = LinearVariationalSolver(problem)
    params = solver.parameters
    params["linear_solver"] = "cg"
    params["preconditioner"] = "hypre_euclid"
    solver.solve()

    # Stabilized Eikonal equation
    print "max cell size:", mesh.hmax()
    eps = Constant(mesh.hmax()/25.)
    F = sqrt(inner(grad(y), grad(y)))*v*dx - f*v*dx + eps*inner(grad(y), grad(v))*dx
    # also works:
    #F = inner(grad(y), grad(y))*v*dx - f*v*dx + eps*inner(grad(y), grad(v))*dx

    problem = NonlinearVariationalProblem(F, y, bc, J=derivative(F, y))
    solver  = NonlinearVariationalSolver(problem)
    params = solver.parameters
    params["newton_solver"]["linear_solver"] = "bicgstab"
    params["newton_solver"]["preconditioner"] = "hypre_euclid"
    solver.solve()
    #solve(F==0, y, bc, solver_parameters=params)
    return y

if __name__ == "__main__":
    y = distance_boundary()
    pughm.Plotter().plot(y, title="distance", interactive=True)
