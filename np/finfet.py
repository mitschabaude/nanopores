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

phys = nanopores.Physics("finfet", geo, xdopants=dopants(Ndop))
#print vars(phys)

pde = nanopores.NonstandardPB(geo, phys)
pde.solve()
pde.visualize()



    

