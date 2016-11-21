# (c) 2016 Gregor Mitscha-Baude
"interpolate 3D diffusivity from piecewise 2D data"
#from nanopores.geometries import pughpore
from nanopores.models import pughpore as pugh
from nanopores.tools.interpolation import harmonic_interpolation
from scipy.interpolate import interp1d
from folders import fields as f
import dolfin
import nanopores
#
#def zsorted(data, field):
#    z = [x[2] for x in data["x"]]
#    J = data[field]
#    I = sorted(range(len(z)), key=lambda k: z[k])
#    z1 = [z[i] for i in I]
#    J1 = [J[i] for i in I]
#    return z1, J1

data = f.get_fields("pugh_diff2D", rMolecule=0.11)
z = [x[2] for x in data["x"]]
data, z = f._sorted(data, z)
Dz = data["D"]

data = f.get_fields("pugh_diff3D_test", rMolecule=0.11, bulkbc=True)
x = [x[0] for x in data["x"]]
data, x = f._sorted(data, x)
D = [d[2][2] for d in data["D"]]

l0 = pugh.pughpore.params["l3"]
r = 0.11
eps = 1e-2
R = l0/2. - r
x += [R, l0/2.+eps]
x = list(reversed([l0/2.-t for t in x]))
x += [100.]

D += [eps, eps]
D00 = D[0]
D = list(reversed([d/D00 for d in D]))
D += [1.]

fz = interp1d(z, Dz)
fx = interp1d(x, D)
import matplotlib.pyplot as plt
plt.plot(z, fz(z), "s-")
plt.figure()
plt.plot(x, fx(x), "s-")

nanopores.add_params(
dim = 3,
h = 1.,
)

functions, mesh = f.get_functions("pugh_distance", h=h)
dist = functions["y"]

#f = lambda x: D0*float(fi(x[dim-1]))
def Diff(x):
    r = dist(x)
    z = x[-1]
    #print "r, z", r, z
    return D0*float(fx(r)*fz(z))

setup = pugh.Setup(dim=dim, x0=None, h=h)
D0 = setup.phys.D

D = harmonic_interpolation(setup, [], [], dict(fluid=Diff), {})
                           #dict(upperb=D0, sideb=D0, lowerb=D0, dnaouterb=D0))
pugh.Plotter(setup).plot(D, title="D")
dolfin.interactive()
plt.show()

"""
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
"""

