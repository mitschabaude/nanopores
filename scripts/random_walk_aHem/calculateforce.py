import math
import numpy as np
from nanopores import *
from nanopores.physics.exittime import ExitTimeProblem
from dolfin import *
def calculateforce(clscale=6., subdomain=None):
    geo_params = dict(
        l3 = 60.,#60
        l4 = 10.,
        R = 60.,#60
        x0 = None, #[5., 0., 10.], # |x0| > 2.2
        exit_i = 1,
    )
    phys_params = dict(
        bV = .5,
        ahemqs = 0.01,
        rTarget = 0.5*nm,
        bulkcon = 1000.,
        T = 273+30
    )
    skip_stokes = True
    StokesProblem.method["iterative"] = True
    taylorhood = True # if True, use P2-P1 discretization for Stokes instead of P1-P1.
    # (True leads too much bigger system but better convergence of iterative solver)
    StokesProblem.method["kparams"].update(
        monitor_convergence = False,
        relative_tolerance = 1e-10,
        absolute_tolerance = 1e-5,
        maximum_iterations = 2000,
        nonzero_initial_guess = True,
    )

    t = Timer("meshing")
    meshdict = generate_mesh(clscale, "aHem", **geo_params)
    print "Mesh generation time:",t.stop()

    t = Timer("reading geometry")
    geo = geo_from_xml("aHem")
    print "Geo generation time:",t.stop()

    phys = Physics("pore_molecule", geo, **phys_params)
        
    pde = PNPS(geo, phys, taylorhood=taylorhood)
    pde.tolnewton = 1e-2
    if skip_stokes:
        pde.solvers.pop("Stokes")
    pde.solve()

    (v, cp, cm, u, p) = pde.solutions(deepcopy=True)
    F, Fel, Fdrag = phys.Forces(v, u)
    
    # save mesh and forces
    File("mesh.xml") << geo.mesh
    File("F.xml") << F
    File("Fel.xml") << Fel
    File("Fdrag.xml") << Fdrag

    for domain in ["pore", "poretop", "porecenter", "porebottom", "fluid_bulk_top", "fluid_bulk_bottom"]:
        print "Average F in %s:"%domain, assemble(F[2]*geo.dx(domain))/assemble(Constant(1.0)*geo.dx(domain))

    return geo.mesh, v
    #VV = VectorFunctionSpace(geo.mesh, "CG", 1)
    #return project(F, VV)
    
def loadforces():
    mesh = Mesh("mesh.xml")
    V = VectorFunctionSpace(mesh, "CG", 1)
    F = Function(V, "F.xml")
    Fel = Function(V, "Fel.xml")
    Fdrag = Function(V, "Fdrag.xml")
    return F, Fel, Fdrag
    
if __name__ == "__main__":
    add_params(scale = 10.)
    mesh, v = calculateforce(clscale=scale)    
