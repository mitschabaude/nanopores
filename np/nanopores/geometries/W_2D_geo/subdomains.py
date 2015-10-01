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
    "ions":"fluid",
    "bulk":{"upperb","lowerb"}, #"rightfluidb"},
    "samb":{"outersamb", "uppersamb", "centersamb", "lowersamb"},
    "aub":"loweraub",
    "sinb":{"lowersinb", "outersinb"},
    "outermembraneb":{"outersamb", "outersinb"},
    "innersamb": {"uppersamb", "centersamb", "lowersamb", },
    "chargedsamb":"innersamb",
    "chargedsinb":{"loweraub", "lowersinb"},
    "membraneb":{"samb", "aub", "sinb"},
    "noslip":{"membraneb","moleculeb"},
    "nopressure":{"upperb","lowerb"},
    "ground":"upperb",
    "bV":"lowerb",
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
    Ry = Rz
    Rx = R
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
                return norm2(x[0], x[1]-x0[2]) <= rMolecule +tolc
            else:
                return False

    # partion membrane into three subdomains
    class MembraneSAM(SubDomain):
        def inside(self, x, on_boundary):
            if lsam < tolc:
                return False
            else:
                return under_line(x[0]-r0, x[1]+l0/2, tan, tolc)  \
                    and between(x[1], (-l0/2 -tolc, l0/2 +tolc))

    class MembraneAu(SubDomain):
        def inside(self, x, on_boundary):
            return under_line(x[0]-rsam, x[1]+l0/2, tan, tolc)  \
                    and between(x[1], ( -l0/2 -tolc, -lsam +l0/2 +tolc))

    class MembraneSiN(SubDomain):
        def inside(self, x, on_boundary):
            return under_line(x[0]-rsin, x[1]+l0/2, tan, tolc)  \
                and between(x[1], ( -l0/2 -tolc, -l0/2 +lsin +tolc))

    # partion pore into three subdomains
    class PoreTop(SubDomain):
        def inside(self, x, on_boundary):
            return over_line(x[0]-r0, x[1]+l0/2, tan, tolc)  \
                and between(x[1], (l0/6 -tolc, l0/2 +tolc))

    class PoreCenter(SubDomain):
        def inside(self, x, on_boundary):
            return over_line(x[0]-r0, x[1]+l0/2, tan, tolc)  \
                and between(x[1], (-l0/6 -tolc, l0/6 +tolc))

    class PoreBottom(SubDomain):
        def inside(self, x, on_boundary):
            return over_line(x[0]-r0, x[1]+l0/2, tan, tolc)  \
                and between(x[1], (-l0/2 -tolc, -l0/6 +tolc))

    return [BulkFluid(), PoreTop(), PoreCenter(), PoreBottom(),
            MembraneSAM(), MembraneAu(), MembraneSiN(), Molecule(),]


def boundaries_list(**params):
    globals().update(params)
    Ry = Rz
    Rx = R

    #define additional variables
    innerfrac = 1 - outerfrac
    sam = None if lsam < tolc or lsam is None else True
    l0 = lsam + lsin + lau

    angle2 = angle/2
    tan = numpy.tan(angle2*numpy.pi/180)
    cos = numpy.cos(angle2*numpy.pi/180)

    r1 = r0 + l0*tan
    rsam = r0 + lsam/cos
    rsin = rsam + rlau

    # exterior fluid boundaries
    class UpperB(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[1], Ry)

    class LowerB(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[1], -Ry)

    class LeftFluidB(SubDomain):
        def inside(self, x, on_boundary):
            if x0 is not None:
                return on_boundary and near(x[0], 0) and  \
                    not between(x[1], (x0[2] - rMolecule +tolc, x0[2] + rMolecule -tolc))
            else:
                return on_boundary and near(x[0], 0)

    class RightFluidB(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], Rx)  \
                and (x[1] < -l0/2 +tolc or x[1] > l0/2 -tolc)

    # Molecule boundaries
    class MoleculeB(SubDomain):
        def inside(self, x, on_boundary):
            if x0 is not None:
                return between(norm2(x[0], x[1]-x0[2]), (rMolecule -tolc, rMolecule +tolc) )
            else:
                return False

    # Membrane boundaries
    class UpperSAMB(SubDomain):
        def inside(self, x, on_boundary):
            return between(x[0], (r1 -tolc, Rx +tolc))  \
                and near(x[1], l0/2)

    class CenterSAMB(SubDomain):
        def inside(self, x, on_boundary):
            return between(tan*(x[1]+l0/2) - (x[0]-r0), (-tolc , +tolc))  \
                and between(x[1], (-l0/2 -tolc, l0/2 +tolc))

    class LowerSAMB(SubDomain):
        def inside(self, x, on_boundary):
            return between(x[0], (r0 -tolc, rsam +tolc)) and near(x[1], -l0/2)

    class LowerAuB(SubDomain):
        def inside(self, x, on_boundary):
            return between(x[0], (rsam -tolc, rsin +tolc)) and near(x[1], -l0/2)

    class LowerSiNB(SubDomain):
        def inside(self, x, on_boundary):
            return between(x[0], (rsin -tolc, Rx +tolc))  \
                and near(x[1], -l0/2)

    class OuterSiNB(SubDomain):
        def inside(self, x, on_boundary):
            return between(x[0], (Rx*innerfrac -tolc, Rx +tolc))  \
                and near(x[1], -l0/2)

    class OuterSAMB(SubDomain):
        def inside(self, x, on_boundary):
            return between(x[0], (Rx*innerfrac -tolc, Rx +tolc))  \
                and near(x[1], l0/2)


    return [UpperB(), LowerB(), LeftFluidB(), RightFluidB(),
            UpperSAMB(), CenterSAMB(), LowerSAMB(), LowerAuB(), LowerSiNB(),
            OuterSiNB(), OuterSAMB(), MoleculeB(),]
