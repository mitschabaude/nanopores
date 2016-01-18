import math
import numpy as np
from nanopores import *
from nanopores.physics.exittime import ExitTimeProblem
from dolfin import *
def calculateforce(clscale=8., tol=1e-2):
    geo_params = dict(
        l3 = 60.,
        l4 = 10.,
        R = 60.,
        x0 = [5., 0., 10.], # |x0| > 2.2
        exit_i = 1,
    )
    phys_params = dict(
        bV = .5,
        ahemqs = 0.01,
        rTarget = 0.5*nm,
        bulkcon = 1000.,
    )

    badexit = {"upperbulkb"}
    goodexit = {"exit"}
    skip_stokes = True


    t = Timer("meshing")
    clscale = 8.
    meshdict = generate_mesh(clscale, "aHem", **geo_params)

    print "Mesh generation time:",t.stop()

    t = Timer("reading geometry")
    geo = geo_from_xml("aHem")

    print "Geo generation time:",t.stop()

    phys = Physics("pore_molecule", geo, **phys_params)

    x0 = geo.params["x0"]
    r0 = math.sqrt(sum(x**2 for x in x0))
    rnear = r0 - geo.params["rMolecule"]
    rfar = r0 + geo.params["rMolecule"]
    xnear = map(lambda x: rnear/r0*x, x0)
    xfar = map(lambda x: rfar/r0*x, x0)

    def avg(u, meas):
        return assemble(u*meas)/assemble(Constant(1.0)*meas)

    def exit_times(tau):
        Tmin = tau(xnear)
        Tmax = tau(xfar)
        Tavg = avg(tau, geo.dS("moleculeb"))
        return (Tmin, Tavg, Tmax)
        
    pde = PNPS(geo, phys)
    pde.tolnewton = 1e-2
    if skip_stokes:
        pde.solvers.pop("Stokes")
    pde.solve()

    (v, cp, cm, u, p) = pde.solutions(deepcopy=True)
    F = phys.Feff(v, u)
    plot(F)
    interactive()
    for domain in ["pore", "poretop", "porecenter", "porebottom", "fluid_bulk_top", "fluid_bulk_bottom"]:
        print "Average F in %s:"%domain, assemble(F[2]*geo.dx(domain))/assemble(Constant(1.0)*geo.dx(domain))

    VV = VectorFunctionSpace(geo.mesh, "CG", 1)
    return project(F, VV)
