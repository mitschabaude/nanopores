" Proof of concept for nanopores.box module "
# FIXME the volume does not get better

from dolfin import *
from nanopores import *
from nanopores.tools.axisym import rotate_z
from nanopores.geometries.curved import Cylinder

# Define cylindrical domain
Rz = 1.
R = 1.
h = 1.
hcross = .2

domain = Box([0., -Rz], [R, Rz])
cross = Box([0., 0.], [R, hcross])
domain.addsubdomains(
    main = domain - cross,
    cross = cross
)
domain.addboundaries(
    lowerb = domain.boundary("bottom"),
    upperb = domain.boundary("top"),
    wall = domain.boundary("right"),
    cross = cross.boundary("bottom"),
)
domain = rotate_z(domain)
geo = domain.create_geometry(lc=h)
# volume of crosssection:
volC = hcross*R**2*pi
areC = R**2*pi
areCh = assemble(Constant(1.)*geo.dS("cross"))
print("Area (approx):", areCh)

# define cylinder
cyl = Cylinder(R=R, L=2.*Rz)
# this causes geo to automatically snap boundaries when adapting
geo.curved = dict(wall = cyl.snap)

plot(geo.submesh("cross"), title=("Mesh 0"))
#plot(geo.mesh, title=("Mesh 0"))
mesh = geo.mesh

num_refinements = 3
for i in range(num_refinements):

    # Mark cells for refinement
    markers = CellFunction("bool", mesh, True)

    # Refine mesh
    mesh = refine(mesh, markers)
    geo.adapt(mesh)
    #geo.snap_to_boundary("wall", cyl.snap)
    """
    cyl.frac = 1 - 0.5*(1. - cyl.frac) # TODO maybe halving is to much

    # Snap boundary
    try:
        mesh.snap_boundary(cyl, True)
    except Exception:
        mesh.snap_boundary(cyl, False)
    """
    areCh = assemble(Constant(1.)*geo.dS("cross"))
    print("Area (approx):", areCh)
    #print "Error A:", abs(areCh - areC)
    
    # Plot mesh
    plot(geo.submesh("cross"), title=("Mesh %d" % (i + 1)))
    #plot(geo.mesh, title=("Mesh %d" % (i + 1)))

interactive()
