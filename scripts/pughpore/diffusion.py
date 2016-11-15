# (c) 2016 Gregor Mitscha-Baude
"investigate ion diffusion constant inside pugh pore by computing drag"
#import numpy as np
import dolfin
import numpy as np
import nanopores as nano
import nanopores.models.pughpore as pugh
import nanopores.physics.simplepnps as pnps
import nanopores.tools.solvers as solvers

default = dict(
    dim = 2,
    h = .6,
    Nmax = 2.7e5,
    rMolecule = 0.152, # radius of K+
    lcMolecule = 0.1,
    H = 200.,
    R = 100.,
    Qmol = 4.,
)

def diffusivity(setup):
    v0 = .001
    geo, phys = setup.geo, setup.phys
    r = setup.geop.rMolecule
    dim = setup.phys.dim
    cyl = setup.phys.cyl

    if geo.mesh.num_cells() < setup.solverp.Nmax:
        pugh.prerefine(setup, True)

    if dim==3:
        pnps.SimpleStokesProblem.method["kparams"]["maximum_iterations"] = 5000
        iterative = True
    else:
        iterative = False
    U0 = dolfin.Constant(tuple(0. for i in range(dim)))
    U1 = dolfin.Constant(tuple((v0 if i==dim-1 else 0.) for i in range(dim)))

    W = pnps.SimpleStokesProblem.space(geo.mesh)
    bcs = [geo.BC(W.sub(0), U0, "dnab"),
           geo.BC(W.sub(0), U0, "memb"),
           #geo.BC(W.sub(0), U0, "sideb"),
           geo.BC(W.sub(0), U1, "moleculeb"),
           geo.BC(W.sub(1), dolfin.Constant(0.0), "upperb")]

    stokes = nano.solve_pde(pnps.SimpleStokesProblem, geo=geo, cyl=cyl,
                            phys=phys, iterative=iterative, bcs=bcs)
    F = stokes.evaluate(phys.Fdrag)["Fdrag"]
    print F

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
        pugh.prerefine(setup, True)

    if dim==3:
        pnps.SimpleStokesProblem.method["kparams"]["maximum_iterations"] = 5000
        iterative = True
    else:
        iterative = False
    U0 = dolfin.Constant(tuple(0. for i in range(dim)))
    W = pnps.SimpleStokesProblem.space(geo.mesh)
    bcs = [geo.BC(W.sub(0), U0, "dnab"),
           geo.BC(W.sub(0), U0, "memb"),
           geo.BC(W.sub(0), U0, "sideb"),
           geo.BC(W.sub(1), dolfin.Constant(0.0), "upperb")]

    gamma = np.zeros((dim, dim))
    for i0 in range(dim):
        U1 = dolfin.Constant(tuple((v0 if i==i0 else 0.) for i in range(dim)))
        bcmol = [geo.BC(W.sub(0), U1, "moleculeb")]
        stokes = nano.solve_pde(pnps.SimpleStokesProblem, geo=geo, cyl=cyl,
                            phys=phys, iterative=iterative, bcs=bcs+bcmol)
        F = stokes.evaluate(phys.Fdrag)["Fdrag"]
        gamma[:,i0] = abs(np.array(F)/v0)

    pi = phys.pi
    eta = phys.eta
    kT = phys.kT

    gamma0 = 6.*pi*eta*r*1e-9
    print "gamma (simulation):", gamma
    print "gamma (stokes law):", gamma0
    print
    D = kT*np.linalg.inv(gamma)
    D0 = kT/gamma0
    print "D (simulation):", D
    print "D (stokes law):", D0
    print "Reducing factor due to confinement:", D/D0
    return D/D0

@solvers.cache_forcefield("pugh_diffusivity", default)
def calculate_diffusivity(X, **params):
    _params = dict(default, **params)
    values = []
    for x0 in X:
        setup = pugh.Setup(x0=x0, **_params)
        D = diffusivity(setup)
        values.append(D)
    return dict(D=values)

@solvers.cache_forcefield("pugh_diffusivity2D", default)
def calculate_diffusivity2D(X, **params):
    _params = dict(default, frac=.5, **params)
    _params["dim"] = 2
    values = []
    for x0 in X:
        setup = pugh.Setup(x0=x0, **_params)
        D = diffusivity(setup)
        values.append(D)
    return dict(D=values)

if __name__ == "__main__":
    print calculate_diffusivity2D([[0.,0.,50.]], cache=False, h=4., Nmax=2e4)
    #print calculate_diffusivity([[0.,0.,0.], [0.,0.,30.]], nproc=2)
