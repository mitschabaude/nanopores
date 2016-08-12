"""
mark domains/boundaries with dolfin MeshFunctions
"""

from dolfin import *
from importlib import import_module
from .params_geo import *
import numpy

synonymes = {
    "pore":{"poretop", "porecenter", "porebottom"},
    "fluid":{"bulkfluid","pore"},
    "sin":"membranesin",
    "au":"membraneau",
    "sam":"membranesam",
    "membrane":{"sin","au","sam"},
    "solid":{"membrane", "molecule"},
    "ions":{"fluid"},
    #"ions":{"fluid","dna"},
    #"ions":{"fluid","dna","membrane"},
    "bulk":{"upperb","lowerb","fluidb"},
    "samb":{"outersamb", "uppersamb", "centersamb", "lowersamb"},
    "aub":"loweraub",
    "sinb":{"lowersinb", "outersinb"},
    "outermembraneb":{"outersamb", "outersinb"},
    "innersamb": {"uppersamb", "centersamb", "lowersamb",},
    "chargedsamb":"innersamb",
    "chargedsinb":{"loweraub", "lowersinb"},  
    "membraneb":{"samb", "aub", "sinb"},
    "noslip":{"membraneb","moleculeb"},
    "nopressure":"bulk",
}

def norm2(x, y):
    return numpy.sqrt(x**2 + y**2)

# lists containing subdomain classes, ordering is important: fluid first, molecule last
def subdomain_list(**params):
    globals().update(params)
    
    def over_line(x, y, tan, tolc):
        return tan*y >= x - tolc

    def under_line(x, y, tan, tolc):
        return tan*y <= x + tolc

    #define additional variables
    sam = None if lsam < tolc or lsam is None else True
    l0 = lsam + lsin + lau

    angle2 = angle/2.0
    tan = numpy.tan(angle2*numpy.pi/180)
    cos = numpy.cos(angle2*numpy.pi/180)

    r1 = r0 + l0*tan
    rsam = r0 + lsam/cos
    rsin = rsam + rlau

    # subdomains
    class BulkFluid(SubDomain):
        def inside(self, x, on_boundary):
            return True  # other domains will overwrite

    class Molecule(SubDomain):
        def inside(self, x, on_boundary):
            if x0 is not None:
                return norm2(norm2(x[0]-x0[0], x[1]-x0[1]), x[2]-x0[2])  \
                    <= rMolecule+tolc
            else:
                return False

    # partion membrane into three subdomains
    class MembraneSAM(SubDomain):
        def inside(self, x, on_boundary):
            if lsam < tolc:
                return False
            else:
                return under_line(norm2(x[0],x[1])-r0, x[2]+l0/2, tan, tolc)  \
                    and between(x[2], (-l0/2 -tolc, l0/2 +tolc))

    class MembraneAu(SubDomain):
        def inside(self, x, on_boundary):
            return under_line(norm2(x[0],x[1])-rsam, x[2]+l0/2, tan, tolc)  \
                and between(x[2], ( -l0/2 -tolc, -lsam +l0/2 +tolc))

    class MembraneSiN(SubDomain):
        def inside(self, x, on_boundary):
            return under_line(norm2(x[0],x[1])-rsin, x[2]+l0/2, tan, tolc)  \
                and between(x[2], ( -l0/2 -tolc, -l0/2 +lsin +tolc))

    # partion pore into three subdomains
    class PoreTop(SubDomain):
        def inside(self, x, on_boundary):
            return over_line(norm2(x[0],x[1])-r0, x[2]+l0/2, tan, tolc)  \
                and between(x[2], (l0/6 -tolc, l0/2 +tolc))

    class PoreCenter(SubDomain):
        def inside(self, x, on_boundary):
            return over_line(norm2(x[0],x[1])-r0, x[2]+l0/2, tan, tolc)  \
                and between(x[2], (-l0/6 -tolc, l0/6 +tolc))

    class PoreBottom(SubDomain):
        def inside(self, x, on_boundary):
            return over_line(norm2(x[0],x[1])-r0, x[2]+l0/2, tan, tolc)  \
                and between(x[2], (-l0/2 -tolc, -l0/6 +tolc))

    return [BulkFluid(), PoreTop(), PoreCenter(), PoreBottom(),
            MembraneSAM(), MembraneAu(), MembraneSiN(), Molecule(),]


def boundaries_list(**params):
    globals().update(params)

    #define additional variables
    sam = None if lsam < tolc or lsam is None else True
    l0 = lsam + lsin + lau
    innerfrac = 1 - outerfrac
    
    angle2 = angle/2.0
    tan = numpy.tan(angle2*numpy.pi/180)
    cos = numpy.cos(angle2*numpy.pi/180)

    r1 = r0 + l0*tan
    rsam = r0 + lsam/cos
    rsin = rsam + rlau

    # exterior fluid boundaries
    class UpperB(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and x[2] >= Rz -tolc
    class LowerB(SubDomain):
        def inside(self, x, on_boundary):
             return on_boundary and x[2] <= -Rz +tolc

    class FluidB(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and norm2(x[0],x[1]) >= R -tolc and  \
                (x[2] < -l0/2 +tolc or x[2] > l0/2 -tolc)

    # Molecule boundaries
    class MoleculeB(SubDomain):
        def inside(self, x, on_boundary):
            if x0 is not None:
                return between(norm2(norm2(x[0]-x0[0],x[1]-x0[1]), x[2]-x0[2]),
                               (rMolecule -tolc, rMolecule +tolc) )
            else:
                return False

    # Membrane boundaries
    class UpperSAMB(SubDomain):
        def inside(self, x, on_boundary):
            return between(norm2(x[0],x[1]), (r1 -tolc, R +tolc))   \
                and near(x[2], l0/2)

    class CenterSAMB(SubDomain):
        def inside(self, x, on_boundary):
            return between(tan*(x[2]+l0/2) - (norm2(x[0],x[1])-r0), (-tolc , +tolc))  \
                and between(x[2], (-l0/2 -tolc, l0/2 +tolc))

    class LowerSAMB(SubDomain):
        def inside(self, x, on_boundary):
            return between(norm2(x[0],x[1]), (r0 -tolc, rsam +tolc))  \
                and near(x[2], -l0/2)

    class LowerAuB(SubDomain):
        def inside(self, x, on_boundary):
            return between(norm2(x[0],x[1]), (rsam -tolc, rsin +tolc))  \
                and near(x[2], -l0/2)

    class LowerSiNB(SubDomain):
        def inside(self, x, on_boundary):
            return between(norm2(x[0],x[1]), (rsin -tolc, R +tolc))  \
                and near(x[2], -l0/2)

    class OuterSiNB(SubDomain):
        def inside(self, x, on_boundary):
            return between(norm2(x[0],x[1]), (R*innerfrac -tolc, R +tolc))  \
                and near(x[2], -l0/2)

    class OuterSAMB(SubDomain):
        def inside(self, x, on_boundary):
            return between(norm2(x[0],x[1]), (R*innerfrac -tolc, R +tolc))  \
                and near(x[2], l0/2)

    return [UpperB(), LowerB(), FluidB(),
            UpperSAMB(), CenterSAMB(), LowerSAMB(), LowerAuB(), LowerSiNB(),
            OuterSiNB(), OuterSAMB(), MoleculeB(),]
