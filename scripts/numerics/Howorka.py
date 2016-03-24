import dolfin
from nanopores import *
from nanopores.geometries.curved import Circle
from mysolve import pbpnps

__all__ = ["setup2D", "solve2D"]

geo_name = "H_geo"
phys_name = "pore_molecule"
nm = import_vars("nanopores.geometries.%s.params_geo" %geo_name)["nm"]

add_params(
h = .5,
z0 = 2.*nm,
bV = -0.1,
Nmax = 1e4,
frac = 0.5,
cheapest = False,
)

def geo_params(z0): 
    x0 = None if z0 is None else [0., 0., z0]
    return dict(
x0 = x0,
rMolecule = 0.5*nm,
moleculeblayer = False,
membraneblayer = False,
)

def phys_params(bV): return dict(
Membraneqs = -0.0,
Qmol = -1.*qq,
bulkcon = 3e2,
dnaqsdamp = .25,
bV = bV,
)

IllposedLinearSolver.stab = 1e0
IllposedNonlinearSolver.newtondamp = 1.
PNPSAxisym.tolnewton = 1e-4
PNPS.tolnewton = 1e-4

def setup2D(**params):
    global z0
    globals().update(params)
    if z0 is not None:
        z0 = round(z0, 4)
    geop = geo_params(z0)
    physp = phys_params(bV)
    
    generate_mesh(h, geo_name, **geop)
    geo = geo_from_name(geo_name, **geop)
    
    if z0 is not None:
        molec = Circle(R=geo.params["rMolecule"], center=geo.params["x0"][::2])
        geo.curved = dict(moleculeb = molec.snap)
    
    phys = Physics(phys_name, geo, **physp)
    return geo, phys

def solve2D(geo, phys, **params):
    globals().update(params)
    return pbpnps(geo, phys, cyl=True, frac=frac, Nmax=Nmax, cheapest=cheapest)

