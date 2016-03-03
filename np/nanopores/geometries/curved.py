""" some parametrized dolfin SubDomains for snapping mesh to curved boundaries """

from dolfin import *

class Cylinder(SubDomain):
    # cylinder in z direction
    # actually the curved part of its boundary

    def __init__(self, R, L, center=(0.,0.,0.), frac=0.85):
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
        if self.inside(x, False):
            x[0] = self.c[0] + (self.R / r)*(x[0] - self.c[0])
            x[1] = self.c[1] + (self.R / r)*(x[1] - self.c[1])
            
            

