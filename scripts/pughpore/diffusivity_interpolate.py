# (c) 2016 Gregor Mitscha-Baude
"interpolate 3D diffusivity from piecewise 2D data"
#from nanopores.geometries import pughpore
from nanopores.models import pughpore as pugh
from nanopores.tools.interpolation import harmonic_interpolation
from scipy.interpolate import interp1d
from folders import fields as f
import dolfin
import nanopores

def zsorted(data, field):
    z = [x[2] for x in data["x"]]
    J = data[field]
    I = sorted(range(len(z)), key=lambda k: z[k])
    z1 = [z[i] for i in I]
    J1 = [J[i] for i in I]
    return z1, J1

#data = f.get_fields("pugh_diffusivity2D", rMolecule=0.152, h=.6, Nmax=2.7e5)
data = f.get_fields("pugh_diffusivity2D", rMolecule=1.)
Z, D = zsorted(data, "D")
fi = interp1d(Z, D)
#import matplotlib.pyplot as plt
#plt.plot(Z, fi(Z))
#plt.show()

nanopores.add_params(
dim = 3,
h = 1.,
)

setup = pugh.Setup(dim=dim, x0=None, h=h)
D0 = setup.phys.D
f = lambda x: D0*float(fi(x[dim-1]))
D = harmonic_interpolation(setup, [], [], dict(pore=f),
                           dict(upperb=D0, sideb=D0, lowerb=D0))#, dnaouterb=D0))
pugh.Plotter(setup).plot(D, title="D")

# solve PNPS and obtain current
geo, phys, solverp = setup.geo, setup.phys, setup.solverp
phys.update(Dp=D, Dm=D)

plotter = pugh.Plotter(setup)
pugh.set_sideBCs(phys, setup.geop, setup.physp)
it = phys.dim==3
pnps = pugh.simplepnps.PNPSFixedPointbV(geo, phys, ipicard=solverp.imax,
           verbose=True, tolnewton=solverp.tol, #taylorhood=True,
           stokesiter=it, iterative=it,
           cyl=phys.cyl) 

print "Number of cells:", geo.mesh.num_cells()
for i in pnps.fixedpoint(ipnp=5):
    v, cp, cm, u, p = pnps.solutions()
    plotter.plot(v, "potential")
    
J = pnps.evaluate(phys.CurrentPNPS)["J"]
print "open pore current:", J, "[A]"

plotter.plot(cp, title="cp")
plotter.plot(cm, title="cm")
dolfin.interactive()


