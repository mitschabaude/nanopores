"""
mark domains/boundaries with dolfin MeshFunctions
"""

from dolfin import *
from math import pow, sqrt
from .params_geo import *

# lists containing subdomain classes, ordering is important fluid first, molecule last
def subdomain_list(**params):
    globals().update(params)
    x0 = params.get('x0')
    rMolecule = params["rMolecule"]
    return [BulkFluid(), DNA(), Membrane(),
            PoreTop(), PoreCenter(), PoreBottom(), Molecule(x0, rMolecule),]
def boundaries_list(**params):
    globals().update(params)
    x0 = params.get('x0')
    rMolecule = params["rMolecule"]
    return [UpperB(), LowerB(), SideB(), ChargedDNAB(), UnchargedDNAB(),
            MembraneB(), CrossTop2D(), CrossCenterTop2D(),
            CrossCenterBottom2D(), CrossBottom2D(), MoleculeB(x0, rMolecule), ]

synonymes = {
    "fluid":{"bulkfluid", "pore"},
    "pore":{"poretop", "porecenter", "porebottom"},
    "ions":"fluid",
    "solid":{"dna", "membrane", "molecule"},
    #"sin":"membrane",
    "lipid":"membrane",
    #"ions":{"fluid", "dna"},
    #"ions":{"fluid", "dna", "membrane"},
    "bulk":{"upperb", "lowerb"}, #, "sideb"},
    #"chargeddnab":{"chargeddnab","unchargeddnab"},
    "dnab":{"chargeddnab", "unchargeddnab"},
    "chargedmembraneb":"membraneb",
    "crosssections2d":{"crosstop2d", "crosscentertop2d", "crosscenterbottom2d", "crossbottom2d"},
    "noslip":{"dnab", "membraneb", "moleculeb"},
    "nopressure":{"upperb","lowerb"},
    "ground":"upperb",
    "bV":"lowerb",
}

def norm2(x, y, z=0.0):
    return sqrt(pow(x,2) + pow(y,2) + pow(z,2))

class BulkFluid(SubDomain):
    def inside(self, x, on_boundary):
        return True  # other domains will overwrite

class Molecule(SubDomain):
    def __init__(self, x0, rMolecule=rMolecule):
        SubDomain.__init__(self)
        self.x0 = x0
        self.rM = rMolecule
    def inside(self, x, on_boundary):
        if self.x0 is not None:
            return norm2(x[0]-self.x0[0], x[1]-self.x0[1], z=x[2]-self.x0[2]) <= self.rM +tolc
        else:
            return False

class DNA(SubDomain):
    def inside(self, x, on_boundary):
        return ( between(norm2(x[0], x[1]), (r0-tolc, r1 +tolc))  \
                 and (abs(x[2])-tolc) <= 0.5*l0 )

class Membrane(SubDomain):
    def inside(self, x, on_boundary):
        return ( between(norm2(x[0], x[1]), (r1 -tolc, R +tolc))  \
                 and abs(x[2]) -tolc <= l1/2 )

# exterior fluid boundaries
class UpperB(SubDomain):
    def inside(self, x, on_boundary):
        return x[2] >= Rz -tolc

class LowerB(SubDomain):
    def inside(self, x, on_boundary):
        return x[2] <= -Rz +tolc

class SideB(SubDomain):
    def inside(self, x, on_boundary):
        return ( norm2(x[0], x[1]) >= R -tolc  \
                 and abs(x[2]) +tolc >= l1/2 )

# Molecule boundaries
class MoleculeB(SubDomain):
    def __init__(self, x0, rMolecule=rMolecule):
        SubDomain.__init__(self)
        self.x0 = x0
        self.rM = rMolecule
    def inside(self, x, on_boundary):
        if self.x0 is not None:
            return between( norm2(x[0]-self.x0[0], x[1]-self.x0[1], z=x[2]-self.x0[2]), (self.rM -tolc, self.rM +tolc) )
        else:
            return False

# DNA boundaries
class ChargedDNAB(SubDomain):
    def inside(self, x, on_boundary):
        #return ( between(norm2(x[0], x[1]), (r0 -tolc, r0 +tolc))  \
        #         or between(norm2(x[0], x[1]), (r1 -tolc, r1 +tolc)) ) \
        #         and between(abs(x[2]), (-tolc, l0/2 +tolc))
        return (( between(norm2(x[0], x[1]), (r0 -tolc, r0 +tolc))  \
                 or between(norm2(x[0], x[1]), (r1 -tolc, r1 +tolc)) )  \
                 and between(abs(x[2]), (l1/2 -tolc, l0/2 +tolc)) ) \
                 or ( between(norm2(x[0], x[1]), (r0 -tolc, r0 +tolc))  \
                 and between(abs(x[2]), (-tolc, l1/2 +tolc)) )

class UnchargedDNAB(SubDomain):
    def inside(self, x, on_boundary):
        return ( between(norm2(x[0], x[1]), (r0 -tolc, r1 +tolc))  \
                 and between(abs(x[2]), (l0/2 -tolc, l0/2 +tolc)) )
        #return ( between(norm2(x[0], x[1]), (r0 -tolc, r0 +tolc))  \
        #         and between(x[2], (-l1/2 -tolc, l1/2 + tolc)) )  \
        #    or ( between(norm2(x[0], x[1]), (r0 -tolc, r1 +tolc))  \
        #         and between(abs(x[2]), (l0/2 -tolc, l0/2 +tolc)) )

# Membrane boundaries
class MembraneB(SubDomain):
    check_midpoint = True
    def inside(self, x, on_boundary):
        return between(norm2(x[0], x[1]), (r1 -tolc, R +tolc)) \
            and near(abs(x[2]), l1/2.) #between(abs(x[2]), (l1/2 -tolc, l1/2 +tolc))

# Center cross-section interfaces (hypersurfaces)
class CrossTop2D(SubDomain):
    def inside(self, x, on_boundary):
        return (near(x[2],l0/2) and between(norm2(x[0],x[1]), (0,r0)))

class CrossCenterTop2D(SubDomain):
    def inside(self, x, on_boundary):
        return (near(x[2],l1/2) and between(norm2(x[0],x[1]), (0,r0)))

class CrossCenterBottom2D(SubDomain):
    def inside(self, x, on_boundary):
        return (near(x[2],-l1/2) and between(norm2(x[0],x[1]), (0,r0)))

class CrossBottom2D(SubDomain):
    def inside(self, x, on_boundary):
        return (near(x[2],-l0/2) and between(norm2(x[0],x[1]), (0,r0)))

# partion pore into three subdomains
class PoreTop(SubDomain):
    def inside(self, x, on_boundary):
        return (between(x[2],(l1/2, l0/2)) and between(norm2(x[0],x[1]), (0,r0)))

class PoreCenter(SubDomain):
    def inside(self, x, on_boundary):
        return (between(x[2],(-l1/2, l1/2)) and between(norm2(x[0],x[1]), (0,r0)))

class PoreBottom(SubDomain):
    def inside(self, x, on_boundary):
        return (between(x[2],(-l0/2,-l1/2)) and between(norm2(x[0],x[1]), (0,r0)))

# get MeshFunctions
def get_subdomain(mesh, x0):
    subdomain = CellFunction("size_t", mesh, 0)
    Fluid().mark(subdomain, 1)
    DNA().mark(subdomain, 2)
    Membrane().mark(subdomain, 3)
    Molecule(x0).mark(subdomain, 4)
    return subdomain

# third argument in mark is for check_midpoint, which defaults to True
def get_boundaries(mesh, x0):
    boundaries = FacetFunction("size_t", mesh, 0)
    UpperB().mark(boundaries, 11, False)
    LowerB().mark(boundaries, 12, False)
    SideB().mark(boundaries, 13, False)
    ChargedDNAB().mark(boundaries, 21, False)
    UnchargedDNAB().mark(boundaries, 22, False)
    MembraneB().mark(boundaries, 31, False)
    MoleculeB(x0).mark(boundaries, 41, False)
    return boundaries

def get_porepartitions(mesh, x0):
    crosssections = get_subdomain(mesh, x0)
    PoreTop().mark(crosssections, 51, False)
    PoreCenter().mark(crosssections, 52, False)
    PoreBottom().mark(crosssections, 53, False)
    Molecule(x0).mark(crosssections, 4)
    return crosssections

def get_crosssections2D(mesh, x0):
    crosssections2D = get_boundaries(mesh, x0)
    CrossTop2D().mark(crosssections2D, 61, False)
    CrossCenterTop2D().mark(crosssections2D, 62, False)
    CrossCenterBottom2D().mark(crosssections2D, 63, False)
    CrossBottom2D().mark(crosssections2D, 64, False)
    MoleculeB(x0).mark(crosssections2D, 41, False)
    return crosssections2D
