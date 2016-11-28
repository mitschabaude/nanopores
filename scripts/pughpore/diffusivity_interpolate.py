# (c) 2016 Gregor Mitscha-Baude
"interpolate 3D diffusivity from piecewise 2D data"
#from nanopores.geometries import pughpore
import numpy as np
from nanopores.models import pughpore as pugh
from nanopores.tools.interpolation import harmonic_interpolation
from eikonal import distance_boundary_from_geo
from scipy.interpolate import interp1d
from folders import fields as f
import dolfin
import nanopores

# read user parameters
params = nanopores.user_params(
    dim = 3,
    h = 6.,
    Nmax = 1e5,
    r = 0.11,
)

# build 1D interpolations from data
data = f.get_fields("pugh_diff2D", rMolecule=params.r)
z = [x[2] for x in data["x"]]
data, z = f._sorted(data, z)
Dz = data["D"]

data = f.get_fields("pugh_diff3D_test", rMolecule=params.r, bulkbc=True)
x = [x[0] for x in data["x"]]
data, x = f._sorted(data, x)
D = [d[2][2] for d in data["D"]]

l0 = pugh.pughpore.params["l3"]
r = params.r
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
#plt.show()

# get geometry, prerefine mesh, compute distance function
setup = pugh.Setup(dim=params.dim, x0=None, h=params.h, Nmax=params.Nmax)
#pugh.prerefine(setup, visualize=True)
dist = distance_boundary_from_geo(setup.geo)
VV = dolfin.VectorFunctionSpace(setup.geo.mesh, "CG", 1)
normal = dolfin.project(dolfin.grad(dist), VV)
plotter = pugh.Plotter(setup)
plotter.plot_vector(normal, "normal")
#functions, mesh = f.get_functions("pugh_distance", h=h)
#dist = functions["y"]

def transformation(n):
    U, S, V = np.linalg.svd(np.matrix(n))
    return V.T

x = [10, 0, 25]
n = normal(x)
print n
T = transformation(n)
print "T inverse", T
print "T", T.T
D = np.diag([1., 2., 3.])
dot = np.dot
print "D transformed", dot(T.T, dot(D, T))
exit()

# interpolate D
def DPore(x):
    r = dist(x)
    z = x[-1]
    return D0*float(fx(r)*fz(z))

def DBulk(x):
    r = dist(x)
    return D0*float(fx(r))

D0 = setup.phys.D
# TODO: save geo along with functions
D = harmonic_interpolation(setup, [], [],
                           dict(bulkfluid=DBulk, poreregion=DPore))
                           #dict(upperb=D0, sideb=D0, lowerb=D0, dnaouterb=D0))
#if not f.exists("Dpugh", **params):
#    f.save_functions("Dpugh", params, dist=dist, D=D)
#    f.update()
plotter.plot(D, title="D")
dolfin.interactive()
exit()
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


