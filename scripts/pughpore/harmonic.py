# (c) 2016 Gregor Mitscha-Baude
"harmonic interpolation -- proof of concept."

from dolfin import *
from nanopores.geometries import pughpore
from nanopores.tools.box import to_mesh

h = 1.
domain = pughpore.get_domain_cyl(h)
domain.compute_entities()
domain.compute_boundaries(True)
domain.entities_to_gmsh(h, True)
domain.physical_to_gmsh(True)


domain.geo = to_mesh()
domain.geo.params = domain.params
domain.geo.import_synonymes(domain.synonymes)

mesh = domain.geo.mesh
#mesh = UnitSquareMesh(10, 10)
V = FunctionSpace(mesh, 'CG', 2)

bdry = CompiledSubDomain('near(x[0], 0)')
value = Expression('x[1]*(1-x[1])', degree=2)
bc = DirichletBC(V, value, bdry, method='pointwise')

# Set quadratic f to quad function
f = Function(V)
f.vector().zero()
bc.apply(f.vector())

# There should be zero error between f and prescribed value on the boundary
facet_f = FacetFunction('size_t', mesh, 0)
bdry.mark(facet_f, 1)
e = (f-value)**2*ds(domain=mesh, subdomain_data=facet_f, subdomain_id=1)

print assemble(e)
plot(f, interactive=True) 