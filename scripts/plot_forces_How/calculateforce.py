import math
import numpy as np
from nanopores import *
#from nanopores.physics.exittime import ExitTimeProblem
from dolfin import *

add_params(
Qmol = -1.,
bV = -0.1,
)

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
        qTarget = lambda Qmol: Qmol,
    )

    t = Timer("meshing")
    generate_mesh(clscale, "H_geo", **geo_params)
    print("Mesh generation time:",t.stop())

    t = Timer("reading geometry")
    geo = geo_from_name("H_geo")
    print("Geo generation time:",t.stop())

    phys = Physics("pore_molecule", geo, **phys_params)
        
    pde = PNPSAxisym(geo, phys)
    pde.tolnewton = 1e-2
    pde.solve()

    (v, cp, cm, u, p) = pde.solutions(deepcopy=True)
    F, Fel, Fdrag = phys.Forces2D(v, u)
    
    # save mesh and forces
    File("mesh2.xml") << geo.mesh
    File("F2.xml") << F
    File("Fel2.xml") << Fel
    File("Fdrag2.xml") << Fdrag

    #return F
    #VV = VectorFunctionSpace(geo.mesh, "CG", 1)
    #return project(F, VV)
    
def loadforces():
    mesh = Mesh("mesh2.xml")
    W = dolfin.FunctionSpace(mesh, "CG", 1)
#    V = dolfin.MixedFunctionSpace((W, W, W))
    V = dolfin.MixedFunctionSpace((W, W))
#    F = Function(V, "F2.xml")
    Fel = Function(V, "Fel2.xml")
    Fdrag = Function(V, "Fdrag2.xml")
    return Fel, Fdrag
    
if __name__ == "__main__":
    add_params(scale = .25)
    calculateforce(clscale=scale)    

