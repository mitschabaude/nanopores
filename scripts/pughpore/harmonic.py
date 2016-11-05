# (c) 2016 Gregor Mitscha-Baude
"harmonic interpolation -- proof of concept."

from dolfin import *
from nanopores.geometries import pughpore
import numpy as np

h = 1.
domain = pughpore.get_domain_cyl(h)
domain.write_gmsh_code(h)

N = 1000
R, H = domain.params["R"], domain.params["H"]
px = R*np.random.random(N)
py = H*np.random.random(N) - H/2.
points = [[x, y, 0.] for x, y in zip(px, py)]
#points = [[(t+1.)/11., (t+1.)/11., 0.] for t in range(200)]

domain.insert_points_2D(points, h, forbidden=["dna", "membrane"])
   
def near_point(x, x0):
    #if all(near(t, t0, 1e-3) for t, t0 in zip(x, x0)):
    #    print "near", x0
    return all(near(t, t0, 1e-3) for t, t0 in zip(x, x0))

class Points(SubDomain):
    def inside(self, x, on_boundary):
        return any(near_point(x, x0) for x0 in points)

geo = domain.code_to_mesh()

geo = domain.geo
mesh = geo.mesh

plot(mesh)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)
a = inner(grad(u), grad(v))*dx
L = Constant(0.)*v*dx(domain=mesh)

# mark vertex function
co = mesh.coordinates()
node_set = set()
for p in points:
    wh = np.where((co[:,0] - p[0])**2 + (co[:,1] - p[1])**2 < 1e-5)[0]
    if wh.shape[0]:
        node_set.add(wh[0])

# Get dofs corresponsing to vertices
dof_set = np.array(vertex_to_dof_map(V)[list(node_set)], dtype='intc')

# Assemble the system
A, b = assemble_system(a, L)

# Boundary value to be prescribed
bc_f = Expression("sin(x[0]/5.)*sin(x[1]/5.)")
# Manual application of bcs
# Modif A: zero bc row & set diagonal to 1
A.ident_local(dof_set)
A.apply('insert')

# Modif b: entry in the bc row is taken from bc_f
bc_values = interpolate(bc_f, V).vector().array()
b_values = b.array()
b_values[dof_set] = bc_values[dof_set]
b.set_local(b_values)
b.apply('insert')

# apply normal bc
bcs = geo.BC(V, Constant(0.), "dnab").bcs + geo.BC(V, Constant(0.), "memb").bcs
[bc.apply(A) for bc in bcs]
[bc.apply(b) for bc in bcs]

# Check that auto matches manual
u = Function(V)
#solve(a==L, u, bcs=bcs)
solve(A, u.vector(), b)
plot(u)
interactive()
