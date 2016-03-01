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

domain = Box([0, -Rz], [Rx, Rz]) #| Box([0+1., -Rz+1.], [Rx+1., Rz+1.])
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
    #fluid = domain - dna
)

domain.addboundaries(
    lowerb = domain.boundary("bottom"),
    upperb = domain.boundary("top"),
    sideb = domain.boundary("right") - membrane.boundary("right"),
    dnab = dna.boundary() - membrane.boundary("left"),
    membraneb = membrane.boundary("top", "bottom"),
)
domain2D = domain
domain2D.plot()
domain = rotate_z(domain, nrot=3)
"""
print "DOMAIN 2D BEFORE"
for item in domain.__dict__.items():
    print "%s : %s" % item
print
print "DOMAIN 3D"
for item in domain.__dict__.items():
    print "%s : %s" % item
print
print "DOMAIN 2D AFTER"
for item in domain2D.__dict__.items():
    print "%s : %s" % item
print
"""
geo = domain.create_geometry(lc=lc)
#plot(geo.mesh)
print "\nBoundaries:\n",geo._physical_boundary
print "\nSubdomains:\n",geo._physical_domain
domain.plot()

for k in 0, 1, 2, 3:
    print "\nEntities[%d]" %k
    for i, (en, st) in enumerate(zip(domain.entities[k], domain.esets[k])):
        print i,":",en,":",st

print "\ndomain:", domain.indexset
if hasattr(domain, "indexsets"):
        print domain.indexsets  
for sub in domain.subdomains + domain.boundaries:
    print sub.name, ":", sub.indexset
    if isinstance(sub, Box):
        print sub.indexsets
        
        

"""


# That was it!


"""
