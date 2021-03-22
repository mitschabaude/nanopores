import dolfin
import matplotlib.pyplot as plt
from nanopores.geometries.pughpore import (params, get_geo_cyl, get_domain,
                                           get_geo)
from nanopores import plot_sliced, user_params
from nanopores.tools.balls import Box, Ball, EmptySet, set_tol, union

up = user_params(params, h=1.)

r = up.rMolecule
lc = up.lcMolecule

positions = [
    [3., 0., 26.],
    [-0.9, .9, 2.],
    [0.1, .8, -24.],
]
molecules = [Ball(x0, r=r, lc=lc) for x0 in positions]

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

solid = membrane | dna | union(molecules)
solid.addsubdomains(membrane=membrane)
solid.addsubdomains(dna=dna)
#solid.addball(molecule, "molecule", "moleculeb")
#solid.addball(molecules[0], "molecule", "moleculeb")
#for i, ball in enumerate(molecules):
#    solid.addball(ball, "molecule%d" %i, "molecule%db" %i)
solid.addballs(molecules, "molecule1", "molecule1b")

solid.addboundaries(
    dnainnerb = domain.getboundary("dnainnerb"),
    dnaupperb = domain.getboundary("dnaupperb"),
    dnaouterb = domain.getboundary("dnaouterb"),
    dnalowerb = domain.getboundary("dnalowerb"),
    memb = domain.getboundary("memb"),
)
print("COMPUTING SOLID")
solidgeo = solid.create_geometry(lc=up.h)
print(solidgeo)

from nanopores.dirnames import DROPBOX
file1 = dolfin.File(DROPBOX + "/geo_pugh.pvd")
file1 << solidgeo.subdomains

plot_sliced(solidgeo, scalarbar=False)
#solidgeo.plot_subdomains()
dolfin.interactive()

#print "COMPUTING DOMAIN"
#geo = get_geo(lc=up.h, **up)
#print geo
#print geo.params
#plot_sliced(geo, scalarbar=False, backend="matplotlib")

