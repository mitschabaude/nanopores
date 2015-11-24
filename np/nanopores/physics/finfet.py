import dolfin
from nanopores.physics.params_physical import *

permSi = eperm*11.7
permOx = eperm*3.9
permDopant = lambda: permSi/4

ni = 1e10/cm**3

Dfin = -1e15/cm**3
Dsource = 7e18/cm**3
Ddopant = lambda: Dsource*5

nDfin = lambda: ni**2/Dfin
nDsource = lambda: ni**2/Dsource
nDdopant = lambda: ni**2/(Dsource*5)

vS = 0.
vD = .5
vG = .1

xdopants = []

synonymes = {
    "si": {"sourcendrain", "fin", "gate"}
}

# use with geo.pwconstBC(V, "bV") --> piece-wise boundary condition
bV = dict(
    sourceb = "vS",
    drainb = "vD",
    gateb = "vG",
)

permittivity = dict(
    oxide = "permOx",
    gate = "permSi",
    fin = "permSi",
    sourcendrain = "permSi",
    dopants = "permDopant",
)

D = dict(
    fin = "Dfin",
    sourcendrain = "Dsource",
    dopants = "Ddopant",
    gate = 0.,
    oxide = 0.,
)

p0 = dict(
    fin = "Dfin",
    sourcendrain = "nDsource",
    dopants = "nDdopant",
    gate = 0.,
    oxide = 0.,
)

n0 = dict(
    fin = "nDfin",
    sourcendrain = "Dsource",
    dopants = "Ddopant",
    gate = 0.,
    oxide = 0.,
)

def add_dopants(geo, xdopants):
    rdop = geo.params["rdop"]
    
    def dist(x, x0):
        return dolfin.sqrt(sum((t-t0)**2 for (t,t0) in zip(x,x0))) 

    class Dopants(dolfin.SubDomain):
        def __init__(self, xi):
            dolfin.SubDomain.__init__(self)
            self.xi = xi
        def inside(self, x, on_boundary):
            return any(dist(x, x0) <= rdop for x0 in self.xi)

    geo.add_subdomain("dopants", Dopants(xdopants))
    return

# TODO: doesn't work    
#def add_synonymes(geo):
#    geo.import_synonymes(synonymes)

def lscale(geo):
    try:
        return geo.parameter("nm")/nm
    except:
        return 1.0
def grad():
    def grad0(u):
        return lscale*dolfin.nabla_grad(u)
    return grad0
def div():
    def div0(u):
        return lscale*dolfin.transpose(dolfin.nabla_div(u))
    return div0
