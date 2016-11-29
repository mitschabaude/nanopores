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
    h = 2.,
    Nmax = 2e6,
    r = 0.11,
)

if not f.exists("Dpugh", **params):
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

    # get geometry, prerefine mesh, compute distance function
    setup = pugh.Setup(dim=params.dim, x0=None, h=params.h, Nmax=params.Nmax,
                cheapest=True)
    pugh.prerefine(setup, visualize=True)
    dist = distance_boundary_from_geo(setup.geo)
    VV = dolfin.VectorFunctionSpace(setup.geo.mesh, "CG", 1)
    normal = dolfin.project(dolfin.grad(dist), VV)
    plotter = pugh.Plotter(setup)
    plotter.plot_vector(normal, "normal")
    #functions, mesh = f.get_functions("pugh_distance", h=h)
    #dist = functions["y"]

    def transformation(n, Dn, Dt):
        D = np.diag([Dn, Dt, Dt]) if len(n)==3 else np.diag([Dn, Dt])
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

    D = lambda i: dict(
        bulkfluid = lambda x: DBulk(x, i),
        poreregion = lambda x: DPore(x, i),
    )

    # TODO: save geo along with functions
    DD = [harmonic_interpolation(setup, subdomains=D(i)) for i in range(params.dim)]
    D = dolfin.Function(VV)
    dolfin.assign(D, DD)

    f.save_functions("Dpugh", params, dist=dist, D=D)
    f.update()
    plt.show()

def get_D(params):
    functions, mesh = f.get_functions("Dpugh", **params)
    D = functions["D"]
    D = dolfin.as_matrix(np.diag([D[i] for i in range(params.dim)]))
    return D

#D = get_D(params)
functions, mesh = f.get_functions("Dpugh", **params)
D = functions["D"]
D = dolfin.as_matrix(np.diag([D[i] for i in range(params.dim)]))
dolfin.plot(functions["D"][0], title="Dx")
dolfin.plot(functions["D"][1], title="Dz")
dolfin.plot(functions["dist"], title="dist")
dolfin.interactive()

exit()
# solve PNPS and obtain current
setup = pugh.Setup(dim=params.dim, x0=None, h=params.h, Nmax=2e5,
                cheapest=True)
if setup.geo.mesh.num_cells() < params.Nmax:
    pugh.prerefine(setup, visualize=True)
geo, phys, solverp = setup.geo, setup.phys, setup.solverp
phys.update(Dp=D, Dm=D)

plotter = pugh.Plotter(setup)
pugh.set_sideBCs(phys, setup.geop, setup.physp)
dim3 = phys.dim==3
large = geo.mesh.num_cells() > 2e5
pnps = pugh.simplepnps.PNPSFixedPointbV(geo, phys, ipicard=solverp.imax,
           verbose=True, tolnewton=solverp.tol, #taylorhood=True,
           stokesiter=dim3 and large, iterative=dim3,
           cyl=phys.cyl)

print "Number of cells:", geo.mesh.num_cells()
for i in pnps.fixedpoint(ipnp=5):
    v, cp, cm, u, p = pnps.solutions()
    #plotter.plot(v, "potential")
    plotter.plot(cm, "cm")

J = pnps.evaluate(phys.CurrentPNPS)["J"]
print "open pore current:", J, "[A]"

plotter.plot(cp, title="cp")
plotter.plot(cm, title="cm")
plotter.plot_vector(u, "u")
dolfin.interactive()

