from nanopores import *
from nanopores.geometries.curved import Circle
from mysolve import pbpnps

__all__ = ["setup2D", "solve2D"]

geo_name = "H_geo"
phys_name = "howorka"
nm = import_vars("nanopores.geometries.%s.params_geo" %geo_name)["nm"]

add_params(
h = .5,
z0 = 2.*nm,
bV = -0.1,
Qmol = -1.,
Nmax = 1e4,
frac = 0.5,
cheapest = False,
dnaqsdamp = .25,
rMolecule = 0.5,
Rx = 12.,
Ry = 12.,
bulkcon = 3e2,
)

def geo_params(z0, rMolecule, Rx, Ry): 
    x0 = None if z0 is None else [0., 0., z0]
    return dict(
x0 = x0,
rMolecule = rMolecule*nm,
moleculeblayer = False,
membraneblayer = False,
Rx = Rx,
Ry = Ry,
)

def phys_params(bV, Qmol, dnaqsdamp, rMolecule, bulkcon): return dict(
Membraneqs = -0.0,
Qmol = Qmol*qq,
bulkcon = bulkcon,
dnaqsdamp = dnaqsdamp,
bV = bV,
qTarget = Qmol*qq,
rTarget = rMolecule*1e-9
)

IllposedLinearSolver.stab = 1e0
IllposedNonlinearSolver.newtondamp = 1.
PNPSAxisym.tolnewton = 1e-4
PNPS.tolnewton = 1e-4

def setup2D(mesh=None, **params):
    # pass mesh argument to prevent mesh regeneration
    global z0
    global geo # to prevent garbage collection of mesh
    globals().update(params)
    if z0 is not None:
        z0 = round(z0, 4)
    geop = geo_params(z0, rMolecule, Rx, Ry)
    physp = phys_params(bV, Qmol, dnaqsdamp, rMolecule, bulkcon)
    
    if mesh is None:
        generate_mesh(h, geo_name, **geop)
        geo = geo_from_name(geo_name, **geop)
    else:
        geo = geo_from_name(geo_name, mesh=mesh, **geop)
    
    if z0 is not None:
        molec = Circle(R=geo.params["rMolecule"], center=geo.params["x0"][::2])
        geo.curved = dict(moleculeb = molec.snap)
    
    phys = Physics(phys_name, geo, **physp)
    #phys.permittivity = {"default": phys.permittivity["water"]}
    return geo, phys

def solve2D(geo, phys, **params):
    globals().update(params)
    return pbpnps(geo, phys, cyl=True, frac=frac, Nmax=Nmax, cheapest=cheapest)    

# evaluate finite-size model for a number of z positions    
def F_explicit(z, **params):
    import dolfin
    values = []
    for z0 in z:
        geo, phys = setup2D(z0=z0, **params)
        dolfin.plot(geo.boundaries, key="b", title="boundaries")
        pb, pnps = solve2D(geo, phys, **params)
        dolfin.plot(geo.boundaries, key="b", title="boundaries")
        values.append(pnps.zforces())
    F, Fel, Fdrag = tuple(zip(*values))
    return F, Fel, Fdrag
     
# evaluate point-size model for a number of z positions
def F_implicit(z, **params):
    geo, phys = setup2D(z0=None, **params)
    pb, pnps = solve2D(geo, phys, **params)
    values = [pnps.zforces_implicit(z0) for z0 in z]
    F, Fel, Fdrag = tuple(zip(*values))
    return F, Fel, Fdrag

# get discrete force fields from point-size model
def F_field_implicit(**params):
    params["z0"] = None
    geo, phys = setup2D(**params)
    pb, pnps = solve2D(geo, phys, **params)
    (v, cp, cm, u, p) = pnps.solutions()
    F, Fel, Fdrag = phys.Forces(v, u)
    return F, Fel, Fdrag