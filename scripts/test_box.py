" Proof of concept for nanopores.box module "
from nanopores import *
from nanopores.geometries.H_geo.params_geo import *

# Redefine 2D Howorka Pore (without molecule)

domain = Box([0, -Rz], [Rx, Rz])
upperhalf = Box([0, 0], [Rx, Rz])

dna = Box([r0, -l0/2], [r1, l0/2])
membrane = Box([r1, -l1/2], [Rx, l1/2])

porebottom = Box([0, -l0/2], [r0, -l1/2])
porecenter = Box([0, -l1/2], [r0, l1/2])
poretop = Box([0, l1/2], [r0, l0/2])
pore = poretop | porecenter | porebottom

bulkfluid = domain - (dna | membrane | pore)
bulkfluid_top = bulkfluid & upperhalf
bulkfluid_bottom = bulkfluid - upperhalf

domain.addsubdomains(
    dna = dna,
    membrane = membrane,
    poretop = poretop,
    porecenter = porecenter,
    porebottom = porebottom,
    bulkfluid_top = bulkfluid_top,
    bulkfluid_bottom = bulkfluid_bottom,
    #bulkfluid = domain - (dna | membrane | pore),
)

domain.addboundaries(
    lowerb = domain.boundary("bottom"),
    upperb = domain.boundary("top"),
    sideb = domain.boundary("right") - membrane.boundary("right"),
    dnab = dna.boundary() - membrane.boundary("left"),
    membraneb = membrane.boundary("top", "bottom"),
    #additionalb = Box(intervals = [(-2*Rx, 2*Rx), (0, 0)]),
)
from dolfin import tic, toc
tic()
geo = domain.create_geometry(lc=.25*lc, mergevols=True)
print "Create geometry:", toc(), "s"
# That was it!

print "\nBoundaries:\n",geo._physical_boundary
print "\nSubdomains:\n",geo._physical_domain
domain.plot()
print
print domain.indexsets
print

for i,en in enumerate(domain.entities[2]):
    print i,":",en
for sub in domain.subdomains:
    print sub.name, ":", sub.indexset
    
print
for i,en in enumerate(domain.entities[1]):
    print i,":",en
for sub in domain.boundaries:
    print sub.name, ":", sub.indexset

