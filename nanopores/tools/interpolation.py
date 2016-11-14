# (c) 2016 Gregor Mitscha-Baude
"harmonic interpolation."
import numpy as np
import dolfin
from nanopores.tools.geometry import PointBC

# TODO: add points in 3D
# TODO: better use iterative solver for laplace

def harmonic_interpolation(setup, points=(), values=(),
                           subdomains=dict(), boundaries=None):
    geo = setup.geo
    phys = setup.phys
    mesh = geo.mesh

    # Laplace equation
    V = dolfin.FunctionSpace(mesh, "CG", 1)
    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)
    a = dolfin.inner(dolfin.grad(u), dolfin.grad(v))*phys.r2pi*geo.dx()
    L = dolfin.Constant(0.)*v*phys.r2pi*geo.dx()

    # Point-wise boundary condition
    bc = PointBC(V, points, values)
    bcs = [bc]

    # Volume boundary conditions
    for sub, f in subdomains.items():
        bc = geo.VolumeBC(V, sub, f)
        bcs.append(bc)

    # Normal boundary conditions
    if boundaries is not None:
        bc = geo.pwBC(V, "", value=boundaries)
        bcs.extend(bc)

    # Assemble, apply bc and solve
    A, b = dolfin.assemble_system(a, L)
    for bc in bcs:
        bc.apply(A)
        bc.apply(b)

    u = dolfin.Function(V)
    if phys.cyl:
        dolfin.solve(A, u.vector(), b, "bicgstab", "hypre_euclid")
    else:
        dolfin.solve(A, u.vector(), b, "cg", "hypre_amg")
    return u

def harmonic_interpolation_simple(mesh, points, values):
    # Laplace equation
    V = dolfin.FunctionSpace(mesh, "CG", 1)
    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)
    a = dolfin.inner(dolfin.grad(u), dolfin.grad(v))*dolfin.dx
    L = dolfin.Constant(0.)*v*dolfin.dx(domain=mesh)

    # Point-wise boundary condition
    bc = PointBC(V, points, values)

    # Assemble, apply bc and solve
    A, b = dolfin.assemble_system(a, L)
    bc.apply(A)
    bc.apply(b)

    u = dolfin.Function(V)
    dolfin.solve(A, u.vector(), b, "cg", "hypre_amg")
    return u

if __name__ == "__main__":
    from nanopores.geometries import pughpore
    from nanopores.tools.utilities import user_params

    h = 1.
    domain = pughpore.get_domain_cyl(h)

    # create points
    N = user_params(N=1000).N
    R, H = pughpore.square2circle(domain.params["R"]), domain.params["H"]
    px = R*np.random.random(N)
    py = H*np.random.random(N) - H/2.
    points = zip(px, py)

    # prepare mesh containing points
    domain.write_gmsh_code(h)
    domain.insert_points(points, h, forbidden=["dna", "membrane"])
    geo = domain.code_to_mesh()
    mesh = geo.mesh

    # interpolate function from point values
    values = np.sin(np.array(points)[:,1]/5.)
    f = lambda x: np.sin(x[1]/5.)
    fexp = dolfin.Expression("sin(x[1]/5.)", domain=mesh, degree=1)

    u = harmonic_interpolation(geo, points, values,
                               dict(bulkfluid_bottom=f),
                               dict(upperb=fexp, sideb=fexp))

    dolfin.plot(mesh)
    dolfin.plot(u)
    dolfin.interactive()
