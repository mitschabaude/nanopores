import dolfin, nanopores
from random import random
from nanopores.geometries.finfet import finfet, dopants

Ndop = 12

t = dolfin.Timer("mesh")
geo = finfet.create_geometry(lc=.8)
print "Mesh generation time:", t.stop()
print "Number of elements:", geo.mesh.num_cells()
print "Number of vertices:", geo.mesh.num_vertices()
#finfet.plot()

phys = nanopores.Physics("finfet", geo, xdopants=dopants(Ndop), vD=.0, vG=.1)
print vars(phys)
#dolfin.plot(geo.subdomains)
#print geo._physical_domain

lpde = nanopores.LinearNonstandardPB(geo, phys)
lpde.solve()
lpde.visualize()
u = lpde.solution

poisson = nanopores.Poisson_(geo, phys)
poisson.solve()
poisson.visualize()

pde = nanopores.NonstandardPB(geo, phys, u=u)
pde.tolnewton = 1e-10
pde.solve()
pde.visualize()


u = pde.solution
ds = geo.dS("crossl")

j0 = u # TODO: change

J0 = dolfin.assemble(j0*ds)
    

