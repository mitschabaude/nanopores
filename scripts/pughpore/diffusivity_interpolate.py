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
Dt = [d[2][2] for d in data["D"]]
Dn = [d[0][0] for d in data["D"]]

l0 = pugh.pughpore.params["l3"]
r = params.r
eps = 1e-2
R = l0/2. - r
x += [R, l0/2.+eps]
x = list(reversed([l0/2.-t for t in x]))
x += [100.]

Dn += [eps, 0.]
Dn = list(reversed([d/Dn[0] for d in Dn]))
Dn += [1.]
Dt += [eps, 0.]
Dt = list(reversed([d/Dt[0] for d in Dt]))
Dt += [1.]

fz = interp1d(z, Dz)
fn = interp1d(x, Dn)
ft = interp1d(x, Dt)
import matplotlib.pyplot as plt
plt.plot(z, fz(z), "s-")
plt.figure()
plt.plot(x, fn(x), "s-")
plt.plot(x, ft(x), "s-")
#plt.show()
#exit()

# get geometry, prerefine mesh, compute distance function
setup = pugh.Setup(dim=params.dim, x0=None, h=params.h, Nmax=params.Nmax)
pugh.prerefine(setup, visualize=True)
dist = distance_boundary_from_geo(setup.geo)
VV = dolfin.VectorFunctionSpace(setup.geo.mesh, "CG", 1)
normal = dolfin.project(dolfin.grad(dist), VV)
plotter = pugh.Plotter(setup)
plotter.plot_vector(normal, "normal")
#functions, mesh = f.get_functions("pugh_distance", h=h)
#dist = functions["y"]

def transformation(n, Dn, Dt):
    D = np.diag([Dn, Dt, Dt])
    U, S, V = np.linalg.svd(np.matrix(n))
    # U, S = 1., |n|
    # V = orthogonal coordinate basis with n/|n| as first vector
    D = np.dot(V, np.dot(D, V.T))
    return np.diag(np.diag(D))

# interpolate D
D0 = setup.phys.D

def DPore(x, i):
    r = dist(x)
    z = x[-1]
    Dz = float(fz(z))
    n = normal(x)
    Dn = D0*float(fn(r))
    Dt = D0*float(ft(r))
    D = transformation(n, Dz*Dn, Dz*Dt)
    return D[i][i]

def DBulk(x, i):
    r = dist(x)
    n = normal(x)
    Dn = D0*float(fn(r))
    Dt = D0*float(ft(r))
    D = transformation(n, Dn, Dt)
    return D[i][i]

D = [dict(
    bulkfluid = lambda x: DBulk(x, i),
    poreregion = lambda x: DPore(x, i),
) for i in 0, 1, 2]

# TODO: save geo along with functions
Dx = harmonic_interpolation(setup, subdomains=D[0])
Dz = harmonic_interpolation(setup, subdomains=D[2])

#if not f.exists("Dpugh", **params):
#    f.save_functions("Dpugh", params, dist=dist, D=D)
#    f.update()
plotter.plot(Dx, title="Dx")
plotter.plot(Dz, title="Dz")
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


