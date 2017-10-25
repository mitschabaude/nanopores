""" some parametrized dolfin SubDomains for snapping mesh to curved boundaries """

from dolfin import *
#import numpy

class Cylinder(SubDomain):
    # cylinder aligned with z axis
    # actually the curved part of its boundary

    def __init__(self, R, L, center=(0.,0.,0.), frac=0.75):
        # R, L ... radius and length
        SubDomain.__init__(self)
        self.R, self.L, self.c, self.frac = R, L, center, frac
        
    def r(self, x):
        return sqrt((x[0] - self.c[0])**2 + (x[1] - self.c[1])**2)
        
    def z(self, x):
        return abs(x[2] - self.c[2])

    def inside(self, x, _):
        return between(self.r(x), (self.frac*self.R, 1./self.frac*self.R)) and (
               between(self.z(x), (0., 0.5*self.L)) )

    def snap(self, x):
        r = self.r(x)
        #if self.inside(x, False):
        x[0] = self.c[0] + (self.R / r)*(x[0] - self.c[0])
        x[1] = self.c[1] + (self.R / r)*(x[1] - self.c[1])
            
#import numpy as np
class Circle(SubDomain):
    # hole

    def __init__(self, R, center=(0.,0.,0.), frac=0.75):
        SubDomain.__init__(self)
        self.R, self.c, self.frac = R, center, frac
        print "DEBUG", self.c, self.R
        
    def r(self, x):
        return sqrt((x[0] - self.c[0])**2 + (x[1] - self.c[-1])**2)

    def inside(self, x, _):
        return between(self.r(x), (0.*self.frac*self.R, 1./self.frac*self.R))
        
    def on_boundary(self, x):
        return self.R - 1e-1 <= self.r(x) <= self.R + 1e-1

    def snap(self, x):
        r = self.r(x)
        #xold = np.copy(x)
        if self.inside(x, False):
            x[0] = self.c[0] + (self.R / r)*(x[0] - self.c[0])
            x[1] = self.c[-1] + (self.R / r)*(x[1] - self.c[-1])
#        
#        rr = np.sqrt(np.sum((xold - x)**2))
#        if rr > 1e-3:
#            print "snap %.4f" % rr, r, xold, "-->", x
            
            
class Sphere(SubDomain):

    def __init__(self, R, center=(0.,0.,0.), frac=0.75):
        # R, L ... radius and length
        SubDomain.__init__(self)
        self.R, self.c, self.frac = R, center, frac
        
    def r(self, x):
        return sqrt(sum((xi - ci)**2 for xi, ci in zip(x, self.c)))

    def inside(self, x, _):
        return between(self.r(x), (self.frac*self.R, 1./self.frac*self.R))
        
    def snap(self, x):
        r = self.r(x)
        x[0] = self.c[0] + (self.R / r)*(x[0] - self.c[0])
        x[1] = self.c[1] + (self.R / r)*(x[1] - self.c[1])
        x[2] = self.c[2] + (self.R / r)*(x[2] - self.c[2])


# associate the actual boundary with every SubDomain that has .on_boundary()            
class Boundary(SubDomain):
    def __init__(self, subdomain):
        SubDomain.__init__(self)
        self.subdomain = subdomain    
    def inside(self, x, _):
        return self.subdomain.on_boundary(x)
                
# snap vertices on a (interior or exterior) boundary to SubDomain
# TODO: this actually needs just the snap function
def snap_to_boundary(geo, name, subdomain):
    mesh = geo.mesh
    # get vertices that lie on this boundary
    # -- mark a CG1 Function with ones on the boundary
    V = FunctionSpace(mesh, 'CG', 1)
    bc = geo.BC(V, 1., name)
    u = Function(V)
    bc.apply(u.vector())

    # Get vertices sitting on boundary
    d2v = dof_to_vertex_map(V)
    vertices_on_boundary = d2v[u.vector() == 1.0]

    # Mark VertexFunction to check
    #
    for v in vertices_on_boundary:
        x = mesh.coordinates()[v]
        subdomain.snap(x)
        mesh.geometry().set(v, x)

    mesh.smooth(1)
    #vf = VertexFunction('size_t', mesh, 0)
    #vf.array()[vertices_on_boundary] = 1.
    #plot(vf)
    #interactive()
            

