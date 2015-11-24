from nanopores.physics.default import *

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

# TODO: geo.pwBC(V, "bV") piece-wise boundary condition
bV = dict(
    sourceb = "vS",
    drainb = "vD",
    gateb = "vG"
)

permittivity = dict(
    oxide = "permSi",
    si = "permSi",
    dopant = "permDopant",
)

D = dict(
    fin = "Dfin",
    sourcendrain = "Dsource",
    dopant = "Ddopant",
    gate = 0.,
    oxide = 0.,
)

p0 = dict(
    fin = "Dfin",
    sourcendrain = "nDsource",
    dopant = "nDdopant",
    gate = 0.,
    oxide = 0.,
)

n0 = dict(
    fin = "nDfin",
    sourcendrain = "Dsource",
    dopant = "Ddopant",
    gate = 0.,
    oxide = 0.,
)
