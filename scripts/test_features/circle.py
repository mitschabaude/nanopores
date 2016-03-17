" Proof of concept for nanopores.box module "

from dolfin import *
from nanopores import *
from nanopores.geometries.curved import Cylinder, Circle, Boundary, snap_to_boundary

# Define rectangular domain
h = .4

domain = Box([0., 0.], [1., 1.])

domain.addsubdomains(
    main = domain
)
domain.addboundaries(
    boundary = domain.boundary(),
)

geo = domain.create_geometry(lc=h)
#boundaries = geo.boundaries
mesh = geo.mesh

# introduce circle
circ = Circle(R=.4, center=(.5,.5), frac=0.9)
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())
subdomains.set_all(0)
circ.mark(subdomains, 1)
plot(subdomains)
submesh = SubMesh(mesh, subdomains, 0)

#plot(mesh, title="initial")
plot(submesh, title="initial")
#plot(boundaries, title="initial")

submesh.snap_boundary(circ, True)
#plot(mesh, title="after snapping")
#plot(submesh, title="after snapping")

# mark circle boundary in geo.boundaries
circb = Boundary(circ)
boundaries = FacetFunction("size_t", submesh, 0)
circb.mark(boundaries, 1)
plot(boundaries, title="after snapping")

subdomains = CellFunction("size_t", submesh, 0)
physb = dict(circle = (1,))
physd = dict(domain = (0,))
geo = Geometry(mesh=submesh, boundaries=boundaries,
    subdomains=subdomains, physical_domain=physd, physical_boundary=physb)
print geo

# Refine mesh
for i in range(4):
    markers = CellFunction("bool", submesh, True)
    submesh = refine(submesh, markers)
    geo.adapt(submesh)
    geo.snap_to_boundary("circle", circ.snap)
    plot(geo.boundaries, title=str(i))

interactive()

