
from dolfin import *
from .params_geo import *
from math import sqrt, pow

def subdomain_list(**params):
    globals().update(params)
    return [Fluid(), Molecule()]
               
def boundaries_list(**params):
    globals().update(params)
    return [Upper(), Lower(), Side(), MoleculeB()]

synonymes = {
    "solid":"molecule",
    "ions":"fluid",
    "noslip":{"upper","lower","side"},
    "nopressure":{"upper"},
    "inflow":{"moleculeb"},
}

synonymes2 = {
    "solid":"molecule",
    "ions":"fluid",
    "noslip":{"side"},
    "nopressure":{"upper"},
    "inflow":{"moleculeb"},
}

synonymes0 = {
    "solid":"molecule",
    "ions":"fluid",
    "inflow":{"upper","lower","side"},
    "nopressure":{"upper"},
    "noslip":{"moleculeb"},
}

def norm2(x, y, z=0.0):
    return sqrt(pow(x,2) + pow(y,2) + pow(z,2))

class Fluid(SubDomain):
    def inside(self, x, on_boundary):
        return True  # other domains will overwrite

class Molecule(SubDomain):
    def inside(self, x, on_boundary):
        return norm2(x[0], x[1]-z0) <= r + tolc

# exterior fluid boundaries
class Upper(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] >= l/2 - tolc

class Lower(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] <= -l/2 + tolc

class Side(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] >= R - tolc

class MoleculeB(SubDomain):
    def inside(self, x, on_boundary):
        return near(norm2(x[0], x[1]-z0), r)
