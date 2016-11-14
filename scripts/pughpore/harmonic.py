# (c) 2016 Gregor Mitscha-Baude
"harmonic interpolation -- proof of concept."
import numpy as np
import dolfin
import nanopores

def harmonic_interpolation(geo, points, values,
                           subdomains=dict(), boundaries=None):
    mesh = geo.mesh

    # Laplace equation
    V = dolfin.FunctionSpace(mesh, "CG", 1)
    u = dolfin.TrialFunction(V)
    v = dolfin.TestFunction(V)
    a = dolfin.inner(dolfin.grad(u), dolfin.grad(v))*dolfin.dx
    L = dolfin.Constant(0.)*v*dolfin.dx(domain=mesh)

    # Point-wise boundary condition
    bc = nanopores.PointBC(V, points, values)
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
    dolfin.solve(A, u.vector(), b)
    return u


#from dolfin import *
import nanopores.geometries.pughpore as pughpore
from nanopores.models.pughpore import Plotter
from nanopores import user_params

h = 4.
domain = pughpore.get_domain(h)

# create random points
N = user_params(N=1000).N
R, H = domain.params["R"], domain.params["H"]
px = R*np.random.random(N)
py = R*np.random.random(N)
pz = H*np.random.random(N) - H/2.
points = zip(px, py, pz)

# prepare mesh containing points
domain.write_gmsh_code(h)
domain.insert_points(points, h, forbidden=["dna", "membrane"])
geo = domain.code_to_mesh()
mesh = geo.mesh

#f = lambda x: np.sin(x[1]/5.)
#V = dolfin.FunctionSpace(mesh, "CG", 1)
#u = dolfin.Function(V)
#bc = VolumeBC(V, geo, "bulkfluid_bottom", f)
#bc.apply(u.vector())


# interpolate function from point values
values = np.sin(np.array(points)[:,2]/5.)
f = lambda x: np.sin(x[2]/5.)
fexp = dolfin.Expression("sin(x[2]/5.)", domain=mesh, degree=1)

u = harmonic_interpolation(geo, points, values,
                           dict(pore=f),
                           dict(upperb=fexp, lowerb=fexp))

dolfin.plot(mesh)
dolfin.plot(u)
Plotter().plot(u)
dolfin.interactive()
