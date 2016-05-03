import math
import numpy as np
from nanopores import *
from nanopores.physics.exittime import ExitTimeProblem
from dolfin import *
def calculateforce(clscale=6., subdomain=None):
    geo_params = dict(
        l3 = 15.,#60
        l4 = 10.,
        R = 15.,#60
        x0 = None, #[5., 0., 10.], # |x0| > 2.2
        exit_i = 1,
    )
    phys_params = dict(
        bV = .5,
        ahemqs = 0.01,
        rTarget = 0.5*nm,
        bulkcon = 1000.,
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
    File("mesh_test.xml") << geo.mesh
    File("F_test.xml") << F
    File("Fel_test.xml") << Fel
    File("Fdrag_test.xml") << Fdrag

    for domain in ["pore", "poretop", "porecenter", "porebottom", "fluid_bulk_top", "fluid_bulk_bottom"]:
        print "Average F in %s:"%domain, assemble(F[2]*geo.dx(domain))/assemble(Constant(1.0)*geo.dx(domain))

    return F
    #VV = VectorFunctionSpace(geo.mesh, "CG", 1)
    #return project(F, VV)
    
def loadforces():
    mesh = Mesh("mesh_test.xml")
    V = VectorFunctionSpace(mesh, "CG", 1)
    F = Function(V, "F_test.xml")
    Fel = Function(V, "Fel_test.xml")
    Fdrag = Function(V, "Fdrag_test.xml")
    return F, Fel, Fdrag
    
if __name__ == "__main__":
    add_params(scale = 10.)
    calculateforce(clscale=scale)    
