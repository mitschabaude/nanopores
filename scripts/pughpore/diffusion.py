# (c) 2016 Gregor Mitscha-Baude
"investigate ion diffusion constant inside pugh pore by computing drag"
#import numpy as np
import dolfin
import nanopores as nano
import nanopores.models.pughpore as pugh
import nanopores.physics.simplepnps as pnps

params = dict(
    rMolecule = 0.152, # radius of K+
    h = 2.,
)
x0 = [0.,0.,-3.]

setup = pugh.Setup(x0=x0, **params)
geo, phys = setup.geo, setup.phys

pnps.SimpleStokesProblem.method["kparams"]["maximum_iterations"] = 5000

W = pnps.SimpleStokesProblem.space(geo.mesh)
bcs = [geo.BC(W.sub(0), dolfin.Constant((0.,0.,0.)), "dnab"),
       geo.BC(W.sub(0), dolfin.Constant((0.,0.,0.)), "memb"),
       geo.BC(W.sub(0), dolfin.Constant((0.,0.,0.)), "sideb"),
       geo.BC(W.sub(0), dolfin.Constant((0.,0.,1.)), "moleculeb"),
       geo.BC(W.sub(1), dolfin.Constant(0.0), "upperb")]
       
stokes = nano.solve_pde(pnps.SimpleStokesProblem, geo=geo, 
                        phys=phys, iterative=True, bcs=bcs)
F = stokes.evaluate(phys.Fdrag)["Fdrag"]
print F
gamma = abs(F[2])
D = phys.kT/gamma
print "D = ", D
u, p = stokes.solutions()
pugh.Plotter(setup).plot_vector(u)
dolfin.interactive()


#params = dict(h=4., Nmax=1e5, rMolecule=)
##params = dict(h=2., Nmax=6e5)
#H = pugh.pughpore.params["H"]
#eps = 5.
#
#ran = np.linspace(-H/2.+eps, H/2.-eps, 12)
#X = [[0.,0.,t] for t in ran]
#
#hpore = pugh.pughpore.params["hpore"]
#ran2 = np.linspace(-hpore/2.-1., -hpore/2.+4., 6)
#X2 = [[0.,0.,t] for t in ran2]
#
#pugh.F_explicit(X+X2, nproc=6, name="pughcenter", **params)
