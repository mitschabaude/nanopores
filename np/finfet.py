import dolfin
from random import random
from nanopores.geometries.finfet import finfet, leftblock, rightblock, rdop

def dist(x, x0):
    return dolfin.sqrt(sum((t-t0)**2 for (t,t0) in zip(x,x0)))

def inside_dopant(x, x0):
    return dist(x, x0) <= rdop + 1e-12

class Dopants(dolfin.SubDomain):
    def __init__(self, xi):
        dolfin.SubDomain.__init__(self)
        self.xi = xi
    def inside(self, x, on_boundary):
        return any(inside_dopant(x, x0) for x0 in self.xi)

def random_dopants(Ndop):
    X = [[random() for i in range(3)] for i in range(Ndop)] # random in [0, 1]**3
    X1 = [[2*x[0] - 1., x[1], x[2]] for x in X if x[0] > .5]
    X2 = [[2*x[0], x[1], x[2]] for x in X if x[0] <= .5]
    def affine(X, box, R):
        return [[r + t*(s-r-R) for t, r, s in zip(x, box.a, box.b)] for x in X]
    X1 = affine(X1, leftblock, rdop)
    X2 = affine(X2, rightblock, rdop)
    return X1 + X2
    
Ndop = 30    

t = dolfin.Timer("mesh")
geo = finfet.create_geometry(lc=.25)
print "Mesh generation time:", t.stop()
print "Number of elements:", geo.mesh.num_cells()
print "Number of vertices:", geo.mesh.num_vertices()
#finfet.plot()

xi = random_dopants(Ndop)
Dopants(xi).mark(geo.subdomains, 1 + max(t[0] for t in geo._physical_domain.values()))
dolfin.plot(geo.submesh("blocks"))
dolfin.interactive()



    

