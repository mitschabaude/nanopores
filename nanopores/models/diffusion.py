# (c) 2016 Gregor Mitscha-Baude
"calculate ion/protein diffusivity constant/tensor by computing drag"
import dolfin
import numpy as np
import nanopores as nano
import nanopores.physics.simplepnps as pnps

def diffusivity(setup, visualize=False):
    v0 = .001
    geo, phys = setup.geo, setup.phys
    r = setup.geop.rMolecule
    dim = setup.phys.dim
    cyl = setup.phys.cyl

    if geo.mesh.num_cells() < setup.solverp.Nmax:
        setup.prerefine(True)

    if dim==3 and geo.mesh.num_cells()>2e5:
        pnps.SimpleStokesProblem.method["kparams"]["maximum_iterations"] = 5000
        iterative = True
    else:
        pnps.SimpleStokesProblem.method["lusolver"] = "default"
        iterative = False

    U0 = dolfin.Constant(tuple(0. for i in range(dim)))
    U1 = dolfin.Constant(tuple((v0 if i==dim-1 else 0.) for i in range(dim)))

    W = pnps.SimpleStokesProblem.space(geo.mesh)
    bcs = [geo.BC(W.sub(0), U0, "dnab"),
           geo.BC(W.sub(0), U0, "memb"),
           geo.BC(W.sub(0), U0, "sideb"),
           geo.BC(W.sub(0), U1, "moleculeb"),
           geo.BC(W.sub(1), dolfin.Constant(0.0), "upperb")]
    if setup.physp.bulkbc:
        bcs.append(geo.BC(W.sub(0), U0, "bulk"))

    stokes = nano.solve_pde(pnps.SimpleStokesProblem, geo=geo, cyl=cyl,
                            phys=phys, iterative=iterative, bcs=bcs)
    F = stokes.evaluate(phys.Fdrag)["Fdrag"]
    print F
    if visualize:
        stokes.visualize("fluid")

    pi = phys.pi
    eta = phys.eta
    kT = phys.kT

    gamma = abs(F[dim-1]/v0)
    gamma0 = 6.*pi*eta*r*1e-9
    print "gamma (simulation):", gamma
    print "gamma (stokes law):", gamma0
    print
    D = kT/gamma
    D0 = kT/gamma0
    print "D (simulation):", D
    print "D (stokes law):", D0
    print "Reducing factor due to confinement:", D/D0
    return D/D0

def diffusivity_tensor(setup):
    v0 = .001
    geo, phys = setup.geo, setup.phys
    r = setup.geop.rMolecule
    dim = setup.phys.dim
    cyl = setup.phys.cyl

    if geo.mesh.num_cells() < setup.solverp.Nmax:
        setup.prerefine(True)

    if dim==3 and geo.mesh.num_cells()>2e5:
        pnps.SimpleStokesProblem.method["kparams"]["maximum_iterations"] = 5000
        iterative = True #False
    else:
        pnps.SimpleStokesProblem.method["lusolver"] = "default"
        iterative = False

    U0 = dolfin.Constant(tuple(0. for i in range(dim)))
    W = pnps.SimpleStokesProblem.space(geo.mesh)
    bcs = [geo.BC(W.sub(0), U0, "dnab"),
           geo.BC(W.sub(0), U0, "memb"),
           geo.BC(W.sub(0), U0, "sideb"),
           geo.BC(W.sub(1), dolfin.Constant(0.0), "upperb")]
    if setup.physp.bulkbc:
        bcs.append(geo.BC(W.sub(0), U0, "bulk"))

    gamma = np.zeros((dim, dim))
    for i0 in range(dim):
        U1 = dolfin.Constant(tuple((v0 if i==i0 else 0.) for i in range(dim)))
        bcmol = [geo.BC(W.sub(0), U1, "moleculeb")]
        stokes = nano.solve_pde(pnps.SimpleStokesProblem, geo=geo, cyl=cyl,
                            phys=phys, iterative=iterative, bcs=bcs+bcmol)
        F = stokes.evaluate(phys.Fdrag)["Fdrag"]
        gamma[:,i0] = abs(np.array(F)/v0)
        #dolfin.plot(stokes.solutions()[0])
        #dolfin.plot(stokes.solutions()[1])

    pi = phys.pi
    eta = phys.eta
    kT = phys.kT

    gamma0 = 6.*pi*eta*r*1e-9
    print "gamma (simulation):\n", gamma
    print "gamma (stokes law):", gamma0
    print
    D = kT*np.linalg.inv(gamma)
    D0 = kT/gamma0
    print "D (simulation):\n", D
    print "D (stokes law):", D0
    print "Reducing factor due to confinement:\n", D/D0
    return D/D0

if __name__ == "__main__":
    #from nanopores.models.pughpore import Setup
    from nanopores.models.Howorka import Setup
    #setup = Setup(dim=2, Nmax=1e4, h=1., x0=[0.,0.,4.6], dnaqsdamp=0.1)
    setup = Setup(dim=3, Nmax=1.7e5, h=.7, x0=[0.,0.,4.6], dnaqsdamp=0.1)
    diffusivity(setup, True)
    dolfin.interactive()
