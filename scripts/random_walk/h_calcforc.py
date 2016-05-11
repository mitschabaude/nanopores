import math
import numpy as np
from nanopores import *
#from nanopores.physics.exittime import ExitTimeProblem
from dolfin import *
h = .5
z0 = 2.*nm
bV = -.1
Qmol = -1.
Nmax = 1e4
frac = 0.5
cheapest = False
def calculateforce(clscale=6., subdomain=None):
    geo_params = dict(
        x0 = None, 
        rMolecule=0.5*nm,
        moleculeblayer=False,
        membraneblayer=False,
    )
    phys_params = dict(
        Membraneqs = -0.0,
        Qmol = Qmol*qq,
        bulkcon = 3e2,
        dnaqsdamp = .25,
        bV = bV,
    )

    t = Timer("meshing")
    meshdict = generate_mesh(clscale, "H_geo", **geo_params)
    print "Mesh generation time:",t.stop()

    t = Timer("reading geometry")
    geo = geo_from_name("H_geo")
    print "Geo generation time:",t.stop()

    phys = Physics("pore_molecule", geo, **phys_params)
        
    pde = PNPS(geo, phys)
    pde.tolnewton = 1e-2
    pde.solve()

    (v, cp, cm, u, p) = pde.solutions(deepcopy=True)
    Der0,Der1 = phys.Forces2D(v, u)
    
    # save mesh and forces
    File("mesh2.xml") << geo.mesh
#    File("F2.xml") << F
    File("Der0.xml") << Der0
    File("Der1.xml") << Der1
#    File("Fdrag2.xml") << Fdrag

    #return F
    #VV = VectorFunctionSpace(geo.mesh, "CG", 1)
    #return project(F, VV)
    
def h_loadforces():
    mesh = Mesh("mesh2.xml")
    V = VectorFunctionSpace(mesh, "CG", 1)
#    F = Function(V, "F2.xml")
    Der0=Function(V, "Der0.xml")
    Der1=Function(V, "Der1.xml")
#    Fdrag = Function(V, "Fdrag2.xml")
    return Der0,Der1
    
if __name__ == "__main__":
    add_params(scale = 10.)
    calculateforce(clscale=scale)    

