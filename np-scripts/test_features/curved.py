" Proof of concept for nanopores.box module "
# TODO make Cylinder.frac dynamic and refine only near the boundary

from dolfin import *
from nanopores import *
from nanopores.tools.axisym import rotate_z
from nanopores.geometries.curved import Cylinder

# Define cylindrical domain
Rz = 1.
R = 1.
h = .5

domain = Box([0., -Rz], [R, Rz])
domain.addsubdomains(
    main = domain
)
domain.addboundaries(
    lowerb = domain.boundary("bottom"),
    upperb = domain.boundary("top"),
    wall = domain.boundary("right"),
)
domain = rotate_z(domain)
geo = domain.create_geometry(lc=h)
mesh = geo.mesh

# add cylinder boundary subdomain
cyl = Cylinder(R=R, L=2.*Rz, frac=0.9)
mesh.snap_boundary(cyl)

# Refine and snap mesh
print "N =", mesh.num_cells()
plot(mesh, title="Mesh 0")

num_refinements = 3
for i in range(num_refinements):

    # Mark cells for refinement
    markers = MeshFunction("bool", mesh, mesh.topology().dim())
    markers.set_all(True)

    # Refine mesh
    mesh = refine(mesh, markers)
    cyl.frac = 1 - 0.5*(1. - cyl.frac)

    # Snap boundary
    try:
        mesh.snap_boundary(cyl, True)
    except Exception:
        mesh.snap_boundary(cyl, False)

    # Plot mesh
    plot(mesh, title=("Mesh %d" % (i + 1)))

interactive()
