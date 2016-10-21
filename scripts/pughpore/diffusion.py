# (c) 2016 Gregor Mitscha-Baude
"investigate ion diffusion constant inside pugh pore by computing drag"
#import numpy as np
import dolfin
import nanopores as nano
import nanopores.models.pughpore as pugh
import nanopores.physics.simplepnps as pnps
import nanopores.tools.solvers as solvers

default = dict(
    h = 4.,
    Nmax = 2e5,
    rMolecule = 0.152, # radius of K+
)

def diffusivity(setup):
    v0 = 1.
    geo, phys = setup.geo, setup.phys
    r = setup.geop.rMolecule
    
    goal = lambda v: phys.Fbare(v, 2)
    pb = pnps.SimpleLinearPBGO(geo, phys, goal=goal)
    for i in pb.adaptive_loop(setup.solverp.Nmax):
        dolfin.plot(geo.submesh("solid"), key="b", title="solid mesh")
    
    pnps.SimpleStokesProblem.method["kparams"]["maximum_iterations"] = 5000
    
    W = pnps.SimpleStokesProblem.space(geo.mesh)
    bcs = [geo.BC(W.sub(0), dolfin.Constant((0.,0.,0.)), "dnab"),
           geo.BC(W.sub(0), dolfin.Constant((0.,0.,0.)), "memb"),
           geo.BC(W.sub(0), dolfin.Constant((0.,0.,0.)), "sideb"),
           geo.BC(W.sub(0), dolfin.Constant((0.,0.,v0)), "moleculeb"),
           geo.BC(W.sub(1), dolfin.Constant(0.0), "upperb")]
           
    stokes = nano.solve_pde(pnps.SimpleStokesProblem, geo=geo, 
                            phys=phys, iterative=True, bcs=bcs)
    F = stokes.evaluate(phys.Fdrag)["Fdrag"]
    print F
    
    pi = phys.pi
    eta = phys.eta
    kT = phys.kT
    
    gamma = abs(F[2]/v0)
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

@solvers.cache_forcefield("pugh_diffusivity", default)
def calculate_diffusivity(X, **params):
    _params = dict(default, **params)
    values = []
    for x0 in X:
        setup = pugh.Setup(x0=x0, **_params)
        D = diffusivity(setup)
        values.append(D)
    return dict(D=[D])
    
if __name__ == "__main__":
    print calculate_diffusivity([[0.,0.,0.], [0.,0.,30.]], nproc=2)
