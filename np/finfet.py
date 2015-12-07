import dolfin, nanopores
from dolfin import *
from random import random
from nanopores.geometries.finfet import finfet, dopants
from nanopores import eperm

Ndop = 12

t = dolfin.Timer("mesh")
geo = finfet.create_geometry(lc=.8)
print "Mesh generation time:", t.stop()
print "Number of elements:", geo.mesh.num_cells()
print "Number of vertices:", geo.mesh.num_vertices()
#finfet.plot()
t = dolfin.Timer("init")
phys = nanopores.Physics("finfet", geo, dopants=dopants(Ndop), vD=.2, vG=.2)
phys.add_dopants
print phys
#dolfin.plot(geo.submesh("sourcendrain"))
#print geo._physical_domain

pde = nanopores.NonstandardPB(geo, phys)
pde.tolnewton = 1e-5
pde.newtondamp = 1.
print "PDE initialization time:", t.stop()
t = dolfin.Timer("solve")
pde.solve()
print "PDE solve time:", t.stop()
u = pde.solution
print phys
dolfin.plot(u, title="potential", interactive=True)
#pde.visualize("sourcendrain")
ds = geo.dS("crossl")

j0 = u # TODO: change

J0 = dolfin.assemble(j0*ds)

