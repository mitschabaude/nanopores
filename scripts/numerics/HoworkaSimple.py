import dolfin, math    
from nanopores import *
from nanopores.geometries.curved import Circle
from nanopores.physics.simplepnps import *

__all__ = ["setup2D", "solve2D"]

geo_name = "H_geo"
phys_name = "pore"
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
imaxh = 10,
tol = 1e-5,
)

def geo_params(z0): 
    x0 = None if z0 is None else [0., 0., z0]
    return dict(
x0 = x0,
rMolecule = 0.5*nm,
moleculeblayer = False,
membraneblayer = False,
#Rx = 20.,
#Ry = 20.,
)

def phys_params(bV, Qmol, dnaqsdamp): return dict(
Membraneqs = -0.0,
Qmol = Qmol*qq,
bulkcon = 3e2,
dnaqsdamp = dnaqsdamp,
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
    physp = phys_params(bV, Qmol, dnaqsdamp)
    
    generate_mesh(h, geo_name, **geop)
    geo = geo_from_name(geo_name, **geop)
    
    if z0 is not None:
        molec = Circle(R=geo.params["rMolecule"], center=geo.params["x0"][::2])
        geo.curved = dict(moleculeb = molec.snap)
    
    phys = Physics(phys_name, geo, **physp)
    #phys.permittivity = {"default": phys.permittivity["water"]}
    return geo, phys
    
def phys2D(geo, **params):
    globals().update(params)
    physp = phys_params(bV, Qmol, dnaqsdamp)
    phys = Physics(phys_name, geo, **physp)
    return phys

def solve2D_fixedpoint(geo, phys, **params):
    globals().update(params)
    pnps = PNPSFixedPoint(geo, phys, cyl=True, ipicard=imax,
                          tolnewton=tol, verbose=True, iterative=False)
    for i in pnps.fixedpoint(ipnp=2):
        dolfin.plot(pnps.functions["poisson"], key="vf")
        #dolfin.plot(pnps.functions["stokes"].sub(0), key="uf")
    return pnps.converged
    
def solve2D_hybrid(geo, phys, **params):
    globals().update(params)
    pnps = PNPSHybrid(geo, phys, cyl=True, inewton=1, ipicard=imax,
                      tolnewton=tol, verbose=True, iterative=False)
    for i in pnps.fixedpoint():
        dolfin.plot(pnps.functions["pnp"].sub(0), key="uh")
    return pnps.converged    

def solve2D_fixedpoint_bVscheme(geo, phys, bVstep=0.025, **params):
    globals().update(params)
    pnps = PNPSFixedPointbV(geo, phys, cyl=True, ipicard=imax,
                          tolnewton=tol, verbose=True, iterative=False)
    
    for i in pnps.fixedpoint(bVstep=bVstep):
        dolfin.plot(pnps.functions["poisson"], key="vfv")
        #dolfin.plot(pnps.functions["stokes"].sub(0), key="ufv")
        v = pnps.functions["poisson"]
        print "v0 =", v([0., -10.])
    return pnps.converged
    
def solve2D_hybrid_PB(geo, phys, **params):
    globals().update(params)
    pb = solve_pde(SimplePBProblem, geo, phys, cyl=True, iterative=False, tolnewton=1e-2)
    pnps = PNPSHybrid(geo, phys, v0=pb.solution, cyl=True, inewton=1, ipicard=imax,
                      tolnewton=tol, verbose=True, iterative=False)
    for i in pnps.fixedpoint():
        dolfin.plot(pnps.functions["pnp"].sub(0), key="uhpb")
    return pnps.converged
    
    
