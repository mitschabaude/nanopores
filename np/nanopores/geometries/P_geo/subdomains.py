"""
mark domains/boundaries with dolfin MeshFunctions
"""

from dolfin import *
from importlib import import_module
from math import sqrt, pow

synonymes = {
    "fluid":{"bulkfluid","pore"},
    "pore": "porecenter",
    "sin":"membrane",
    "molecule":"dna",
    "moleculeb":"dnab",
    "solid":{"dna", "membrane"},
    "ions":{"fluid"},
    #"ions":{"fluid","dna"},
    #"ions":{"fluid","dna","membrane"},
    "bulk":{"upperb","lowerb",}, #"rightfluidb"},
    "dnab":{"chargeddnab","unchargeddnab"},
    "chargedmembraneb":{"innermembraneb", "openingmembraneb"},
    "membraneb":{"chargedmembraneb", "outermembraneb"},
    "noslip":{"dnab","membraneb"},
    "nopressure":{"upperb","lowerb"},
    "ground": "upperb",
    "bV": "lowerb",
}

# lists containing subdomain classes, ordering is important
def subdomain_list(**params):
    tmp = vars(import_module('.params_geo','nanopores.geometries.P_geo'))
    params.update({key:tmp[key] for key in tmp  \
                   if (not key in params and not key.startswith("__"))})
    tolc = params["tolc"]

    # subdomains
    class BulkFluid(SubDomain):
        def inside(self, x, on_boundary):
            return True  # other domains will overwrite

    class DNA(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] <= params["DNAradius"] +tolc and abs(x[1]) <= params["DNAlength"]/2 + tolc

    class Membrane(SubDomain):
        def inside(self, x, on_boundary):
            return x[0] >= params["r0"] -tolc and abs(x[1]) <= params["l1"]/2 +tolc \
                and abs(x[1]) -tolc <= x[0]*(params["l1"]-params["l0"])/2/(params["r1"]-params["r0"])+(params["r1"]*params["l0"]-params["r0"]*params["l1"])/2/(params["r1"]-params["r0"])

    class PoreCenter(SubDomain):
        def inside(self, x, on_boundary):
            return (between(x[1],(-params["l0"]/2, params["l0"]/2)) and between(x[0], (params["DNAradius"], params["r0"])))

    return [BulkFluid(), DNA(), Membrane(), PoreCenter(), ]


def boundaries_list(**params):
    tmp = vars(import_module('.params_geo','nanopores.geometries.P_geo'))
    params.update({key:tmp[key] for key in tmp  \
                   if (not key in params and not key.startswith("__"))})
    tolc = params["tolc"]

    # exterior fluid boundaries
    class UpperB(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[1], params["Ry"])

    class LowerB(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[1], -params["Ry"])

    class LeftFluidB(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], 0) and  \
                abs(x[1]) > params["DNAlength"]/2

    class RightFluidB(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], params["Rx"]) and  \
                abs(x[1]) >= params["l1"]/2 -tolc

    # DNA boundaries
    class ChargedDNAB(SubDomain):
        def inside(self, x, on_boundary):
            return between(x[0], (params["DNAradius"] -tolc, params["DNAradius"] +tolc)) and  \
                between(x[1], (-params["DNAlength"]/2 -tolc, params["DNAlength"]/2 +tolc))

    class UnchargedDNAB(SubDomain):
        def inside(self, x, on_boundary):
            return between(x[0], (-tolc, params["DNAradius"] +tolc)) and  \
                between(abs(x[1]), (params["DNAlength"]/2 -tolc, params["DNAlength"]/2 +tolc))

    # Membrane boundaries
    class InnerMembraneB(SubDomain):
        def inside(self, x, on_boundary):
            return between(x[0], (params["r0"] -tolc, params["r0"] +tolc)) and between(x[1], (-params["l0"]/2 -tolc, params["l0"]/2 +tolc))

    class OpeningMembraneB(SubDomain):
        def inside(self, x, on_boundary):
            return between(abs(x[1]),(x[0]*(params["l1"]-params["l0"])/2/(params["r1"]-params["r0"])+(params["r1"]*params["l0"]-params["r0"]*params["l1"])/2/(params["r1"]-params["r0"]) -tolc, x[0]*(params["l1"]-params["l0"])/2/(params["r1"]-params["r0"])+(params["r1"]*params["l0"]-params["r0"]*params["l1"])/2/(params["r1"]-params["r0"]) +tolc))

    class OuterMembraneB(SubDomain):
        def inside(self, x, on_boundary):
            return between(x[0], (params["r1"] -tolc, params["Rx"] +tolc)) and between(abs(x[1]), (params["l1"]/2 -tolc, params["l1"]/2 +tolc))

    # cross-section interfaces
    class CrossCenterTop2D(SubDomain):
        def inside(self, x, on_boundary):
            return (near(x[1],params["l0"]/2) and  \
                    between(x[0],(params["DNAradius"] -tolc, params["r0"] +tolc)))

    class CrossCenterBottom2D(SubDomain):
        def inside(self, x, on_boundary):
            return (near(x[1],-params["l0"]/2) and  \
                    between(x[0],(params["DNAradius"] -tolc, params["r0"] +tolc)))

    return [UpperB(), LowerB(), LeftFluidB(), RightFluidB(),
            ChargedDNAB(), UnchargedDNAB(),
            InnerMembraneB(), OpeningMembraneB(), OuterMembraneB(),
            CrossCenterTop2D(), CrossCenterBottom2D(), ]
