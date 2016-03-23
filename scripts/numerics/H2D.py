from nanopores import *
from nanopores.geometries.curved import Circle

__all__ = ["setup", "solve_pnps"]

geo_name = "H_geo"
phys_name = "pore_molecule"
nm = import_vars("nanopores.geometries.%s.params_geo" %geo_name)["nm"]

add_params(
h = 1.,
z0 = 2.*nm,
bV = -0.1,
Nmax = 1e4,
)

def geo_params(z0): return dict(
x0 = [0., 0., z0],
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

def setup(**params):
    # these lines can surely be turned into some more general magic
    globals().update(params)
    geo_params = geo_params(z0)
    phys_params = phys_params(bV)
    
    generate_mesh(h, geo_name, **geo_params)
    geo = geo_from_name(geo_name, **geo_params)

    molec = Circle(R=geo.params["rMolecule"], center=geo.params["x0"][::2])
    geo.curved = dict(moleculeb = molec.snap)

    phys = Physics(phys_name, geo, **phys_params)
    return geo, phys

def solve_pnps(geo, phys):
    pass
