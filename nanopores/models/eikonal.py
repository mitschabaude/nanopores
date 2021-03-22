"solve eikonal equation to get distance to boundary"
from dolfin import *

def distance_boundary_from_geo_OLD(geo, boundary="dnab"):
    mesh = geo.mesh
    V = FunctionSpace(mesh, "CG", 1)
    v = TestFunction(V)
    u = TrialFunction(V)
    f = Constant(1.0)
    y = Function(V)

    bc = geo.BC(V, Constant(0.), boundary).bcs# + \
         #geo.BC(V, Constant(0.), "memb").bcs

    # Initialization problem to get good initial guess for nonlinear problem:
    F1 = inner(grad(u), grad(v))*dx - f*v*dx
    problem = LinearVariationalProblem(lhs(F1), rhs(F1), y, bc)
    solver = LinearVariationalSolver(problem)
    params = solver.parameters
    params["linear_solver"] = "cg"
    params["preconditioner"] = "hypre_euclid"
    solver.solve()

    # Stabilized Eikonal equation
    print("max cell size:", mesh.hmax())
    eps = Constant(mesh.hmax()/25.)
    #F = sqrt(inner(grad(y), grad(y)))*v*dx - f*v*dx + eps*inner(grad(y), grad(v))*dx
    # also works:
    F = inner(grad(y), grad(y))*v*dx - f*v*dx + eps*inner(grad(y), grad(v))*dx

    problem = NonlinearVariationalProblem(F, y, bc, J=derivative(F, y))
    solver  = NonlinearVariationalSolver(problem)
    params = solver.parameters
    #params["newton_solver"]["linear_solver"] = "bicgstab"
    #params["newton_solver"]["preconditioner"] = "hypre_euclid"
    solver.solve()
    #solve(F==0, y, bc, solver_parameters=params)
    return y

def distance_boundary_from_geo(geo, boundary="dnab",
                               stab=2., eps0=1., epsfact=0.5):
    mesh = geo.mesh
    V = FunctionSpace(mesh, "CG", 1)
    v = TestFunction(V)
    du = TrialFunction(V)
    f = Constant(1.0)

    bc = geo.BC(V, Constant(0.), boundary).bcs + \
         geo.BC(V, Constant(0.), "memb").bcs
    eps = eps0
    epsi = Constant(eps)
    h = CellSize(mesh)
    u = Function(V)

    a = epsi*h*inner(grad(du), grad(v))*dx + 2.*inner(grad(u), grad(du))*v*dx
    L = -epsi*h*inner(grad(u), grad(v))*dx + (f - inner(grad(u), grad(u)))*v*dx

    # add SUPG stabilization terms
    tau = Constant(stab)
    a = a + tau*h*inner(grad(u), grad(du))*inner(grad(u), grad(v))*dx
    L = L + tau*0.5*h*(f - inner(grad(u), grad(u)))*inner(grad(u), grad(v))*dx

    du = Function(V)
    problem = LinearVariationalProblem(a, L, du, bc)
    solver = LinearVariationalSolver(problem)
    if geo.dim() == 3:
        params = solver.parameters
        params["linear_solver"] = "bicgstab"
        params["preconditioner"] = "hypre_euclid"

    while eps > 1e-4:
        print("eps", eps, " ", end=' ')
        solver.solve()
        u.assign(u + du)
        eps *= epsfact
        epsi.assign(eps)

    return u


if __name__ == "__main__":
    import nanopores.geometries.pughpore as pughpore
    import nanopores.models.pughpore as pugh
    from nanopores import user_param
    setup = pugh.Setup(dim=3, h=2., x0=None, Nmax=6e5, cheapest=True)
    setup.prerefine()
    geo = setup.geo
    #geo = pughpore.get_geo_cyl(lc=1., x0=None)
    y = distance_boundary_from_geo(geo)
    setup.plot(y, title="distance", interactive=True)
    from numpy import linspace
    from matplotlib.pyplot import plot, show
    t = linspace(0., 20., 100)
    def point(t):
        x = [0.]*setup.geop.dim
        x[0] = t
        return x

    plot(t, [y(point(t0)) for t0 in t], ".-")
    show()
