" Proof of concept for nanopores.box module "
from nanopores import *
from nanopores.tools.axisym import *
from dolfin import plot

# Redefine 2D Howorka Pore (without molecule)
Rz = 10.
Rx = 8.
l0 = 15.
l1 = 2.
r0 = 1.
r1 = 2.5
lc = 1.

domain = Box([0, -Rz], [Rx, Rz])
dna = Box([r0, -l0/2], [r1, l0/2])
membrane = Box([r1, -l1/2], [Rx, l1/2])

porebottom = Box([0, -l0/2], [r0, -l1/2])
porecenter = Box([0, -l1/2], [r0, l1/2])
poretop = Box([0, l1/2], [r0, l0/2])
pore = poretop | porecenter | porebottom

domain.addsubdomains(
    dna = dna,
    membrane = membrane,
    poretop = poretop,
    porecenter = porecenter,
    porebottom = porebottom,
    fluid = domain - (dna | membrane | pore),
)

domain.addboundaries(
    lowerb = domain.boundary("bottom"),
    upperb = domain.boundary("top"),
    sideb = domain.boundary("right") - membrane.boundary("right"),
    dnab = dna.boundary() - membrane.boundary("left"),
    membraneb = membrane.boundary("top", "bottom"),
)

domain = rotate_z(domain)
geo = domain.create_geometry(lc=lc)
#plot(geo.mesh)
print "\nBoundaries:\n",geo._physical_boundary
print "\nSubdomains:\n",geo._physical_domain
domain.plot()

for k in 0, 1, 2, 3:
    print "Entities[%d]" %k
    for i, (en, st) in enumerate(zip(domain.entities[k], domain.esets[k])):
        print i,":",en,":",st
    print

print "domain:", domain.indexset
if hasattr(domain, "indexsets"):
        print domain.indexsets  
for sub in domain.subdomains + domain.boundaries:
    print sub.name, ":", sub.indexset
    if isinstance(sub, Box):
        print sub.indexsets
        
        

"""


# That was it!


"""
