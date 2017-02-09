import dolfin
import matplotlib.pyplot as plt
from nanopores.geometries.pughpore import (params, get_geo_cyl, get_domain,
                                           get_geo)
from nanopores import plot_sliced, user_params
up = user_params(params, h=4.)

#geo2D = get_geo_cyl(lc=1., **up)
#dolfin.plot(geo2D.subdomains, title="subdomains", backend="matplotlib")
#dolfin.plot(geo2D.boundaries, title="boundaries", backend="matplotlib")
#dolfin.plot(geo2D.mesh, backend="matplotlib")
#print geo2D
#dolfin.interactive()

domain = get_domain(**up)
membrane = domain.getsubdomain("membrane")
dna = domain.getsubdomain("dna")
molecule = domain.getsubdomain("molecule")

solid = membrane | dna | molecule
solid.addsubdomains(membrane=membrane)
solid.addsubdomains(dna=dna)
solid.addball(molecule, "molecule", "moleculeb")
solid.addboundaries(
    dnainnerb = domain.getboundary("dnainnerb"),
    dnaupperb = domain.getboundary("dnaupperb"),
    dnaouterb = domain.getboundary("dnaouterb"),
    dnalowerb = domain.getboundary("dnalowerb"),
    memb = domain.getboundary("memb"),
)
print "COMPUTING SOLID"
solidgeo = solid.create_geometry(lc=up.h)
print solidgeo
#dolfin.plot(solidgeo.boundaries, title="boundaries", scalarbar=False)
plot_sliced(solidgeo, scalarbar=False)

#print "COMPUTING DOMAIN"
#geo = get_geo(lc=up.h, **up)
#print geo
#print geo.params
#plot_sliced(geo, scalarbar=False, backend="matplotlib")

dolfin.interactive()