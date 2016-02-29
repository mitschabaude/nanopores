import dolfin, nanopores
from dolfin import *
from random import random
from nanopores.geometries.finfet import finfet, dopants
from nanopores import eperm, cm, qq

Ndop = 12

# --- create mesh and geometrical/physical context information
t = dolfin.Timer("mesh")
geo = finfet.create_geometry(lc=2.)
print "Mesh generation time:", t.stop()
print "Number of elements:", geo.mesh.num_cells()
print "Number of vertices:", geo.mesh.num_vertices()
#finfet.plot()
t = dolfin.Timer("init")
phys = nanopores.Physics("finfet", geo, dopants=dopants(Ndop), vD=.2, vG=.2)
phys.add_dopants
#print phys
#dolfin.plot(geo.submesh("sourcendrain"))
#print geo._physical_domain

# --- definition and solution of PDE
pde = nanopores.NonstandardPB(geo, phys)
pde.tolnewton = 1e-5
pde.newtondamp = 1.
print "PDE initialization time:", t.stop()
t = dolfin.Timer("solve")
pde.solve()
print "PDE solve time:", t.stop()
u = pde.solution

# --- calculation of the current
p0 = geo.pwconst("p0")
n0 = geo.pwconst("n0")
beta = 1./phys.UT
grad = phys.grad

n = n0*exp(beta*u)
p = p0*exp(-beta*u)
mun = 1000*cm**2
mup = 100*cm**2
L = geo.params["length"]
Lfin = geo.params["lfin"]
E = (phys.vD - phys.vS)/L

jn = Constant(mun*qq*E)*n
jp = Constant(mup*qq*E)*p
dx = geo.dx("fin")

Jn = dolfin.assemble(jn*Constant(1./Lfin)*dx)
Jp = dolfin.assemble(jp*Constant(1./Lfin)*dx)

print "Jn [A]: ",Jn
print "Jp [A]: ",Jp

# --- visualization
dolfin.plot(u, title="potential")
dolfin.interactive()

