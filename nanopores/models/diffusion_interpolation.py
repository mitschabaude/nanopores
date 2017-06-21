# (c) 2016 Gregor Mitscha-Baude
"""interpolate 2D/3D diffusivity fields from piecewise 1D profiles.
involves solving the eikonal equation to obtain the distance to boundary."""

import numpy as np
from scipy.interpolate import interp1d
import dolfin

import nanopores
from nanopores.tools.interpolation import harmonic_interpolation
from nanopores.models.eikonal import distance_boundary_from_geo
from nanopores.tools import fields

def transformation(n, Dn, Dt):
    # TODO RR =
    D = np.diag([Dn, Dt, Dt]) if len(n)==3 else np.diag([Dn, Dt])
    U, S, V = np.linalg.svd(np.matrix(n))
    # U, S = 1., |n|
    # V = orthogonal coordinate basis with n/|n| as first vector
    D = np.dot(V, np.dot(D, V.T))
    return np.diag(np.diag(D))

def preprocess_Dr(data, r, normalize=True):
    x = data["x"] #[xx[0] for xx in data["x"]]
    data, x = fields._sorted(data, x)
    Dt = [d[2][2] for d in data["D"]]
    Dn = [d[0][0] for d in data["D"]]

    eps = 1e-2
    x = [-1., -eps, r] + x + [100.]
    if normalize:
        Dn = [d/Dn[-1] for d in Dn]
        Dt = [d/Dt[-1] for d in Dt]

    Dn = [0., 0., eps] + Dn + [1.]
    Dt = [0., 0., eps] + Dt + [1.]

    fn = interp1d(x, Dn)
    ft = interp1d(x, Dt)
    return fn, ft

def Dt_plane(h, r):
    x = r/h
    return 1. - 9./16.*x + 1./8.*x**3 - 45./256.*x**4 - 1/16.*x**5

sinh = np.sinh
acosh = np.arccosh
def Dn_plane(l, r, N=100):
    alpha = acosh(l/r)
    s = 0.
    for n in range(1, N):
        n = float(n)
        K = n*(n+1)/(2*n-1)/(2*n+3)
        s += K*((2*sinh((2*n+1)*alpha)+(2*n+1)*sinh(2*alpha))/(4*(sinh((n+.5)*alpha))**2-(2*n+1)**2*(sinh(alpha))**2) - 1)
    return 1./((4./3.)*sinh(alpha)*s)

def matrix(d):
    return [[d[0], 0., 0.], [0., d[1], 0.], [0., 0., d[2]]]

def diff_profile_plane(r):
    eps = 1e-8
    x = np.linspace(r+eps, r*16., 100)

    Dn = np.array([Dn_plane(xx, r) for xx in x])
    Dt = Dt_plane(x, r)

    data = dict(x=list(x), D=map(matrix, zip(Dn, Dt, Dt)))
    return data

def diff_profile_trivial(r):
    eps = 1e-8
    x = np.linspace(r+eps, r*16., 100)
    D = [1.]*len(x)
    data = dict(x=list(x), D=map(matrix, zip(D, D, D)))
    return data

def diffusivity_field(setup, r, ddata_r=None, ddata_z=None,
                      boundary="dnab", poreregion="poreregion"):
    "interpolates diffusivity field defined on the geometry given by setup"
    "needs: molecule radius r, radial and vertical D profiles ddata_r, ddata_z"
    "the geometry needs to have subdomains poreregion and bulkfluid"

    if ddata_r is None:
        ddata_r = diff_profile_plane(r)

    if ddata_z is None:
        return diffusivity_field_simple(setup, r, ddata_r, boundary=boundary)

    # build 1D interpolations from data
    data = ddata_z
    z = [x[2] for x in data["x"]]
    data, z = fields._sorted(data, z)
    Dz = data["D"]
    fz = interp1d(z, Dz)

    fn, ft = preprocess_Dr(ddata_r, r)

    setup.prerefine(visualize=True)
    dist = distance_boundary_from_geo(setup.geo, boundary)
    VV = dolfin.VectorFunctionSpace(setup.geo.mesh, "CG", 1)
    normal = dolfin.project(dolfin.grad(dist), VV)

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

    D = lambda i: {
        "bulkfluid": lambda x: DBulk(x, i),
        poreregion: lambda x: DPore(x, i)}

    dim = setup.geop.dim
    DD = [harmonic_interpolation(setup, subdomains=D(i)) for i in range(dim)]
    D = dolfin.Function(VV)
    dolfin.assign(D, DD)

    return dict(dist=dist, D=D)

def diffusivity_field_alt(setup, r, ddata_pore, ddata_bulk):
    "interpolates diffusivity field defined on the geometry given by setup"

    # build 1D interpolations from data
    fn, ft = preprocess_Dr(ddata_bulk, r)
    fn_pore, ft_pore = preprocess_Dr(ddata_pore, r, normalize=False)

    setup.prerefine(visualize=True)
    dist = distance_boundary_from_geo(setup.geo)
    VV = dolfin.VectorFunctionSpace(setup.geo.mesh, "CG", 1)
    normal = dolfin.project(dolfin.grad(dist), VV)

    phys = setup.phys
    D0 = phys.kT / (6.*np.pi*phys.eta*r*1e-9)

    def DPore(x, i):
        r = dist(x)
        n = normal(x)
        Dn = D0*float(fn_pore(r))
        Dt = D0*float(ft_pore(r))
        D = transformation(n, Dn, Dt)
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
        nearpore = lambda x: DBulk(x, i),
        poreenter = lambda x: DBulk(x, i),
        porecurrent = lambda x: DPore(x, i),
        porerest = lambda x: DPore(x, i),
    )

    dim = setup.geop.dim
    DD = [harmonic_interpolation(setup, subdomains=D(i)) for i in range(dim)]
    D = dolfin.Function(VV)
    dolfin.assign(D, DD)

    return dict(dist=dist, D=D)

def diffusivity_field_simple(setup, r, ddata_bulk, boundary="poresolidb"):
    "interpolates diffusivity field defined on the geometry given by setup"

    # build 1D interpolations from data
    fn, ft = preprocess_Dr(ddata_bulk, r)

    setup.prerefine(visualize=True)
    dist = distance_boundary_from_geo(setup.geo, boundary)
    VV = dolfin.VectorFunctionSpace(setup.geo.mesh, "CG", 1)
    normal = dolfin.project(dolfin.grad(dist), VV)

    D0 = setup.phys.D

    def DBulk(x, i):
        r = dist(x)
        n = normal(x)
        Dn = D0*float(fn(r))
        Dt = D0*float(ft(r))
        D = transformation(n, Dn, Dt)
        return D[i][i]

    D = lambda i: dict(
        fluid = lambda x: DBulk(x, i),
    )

    dim = setup.geop.dim
    DD = [harmonic_interpolation(setup, subdomains=D(i)) for i in range(dim)]
    D = dolfin.Function(VV)
    dolfin.assign(D, DD)

    return dict(dist=dist, D=D)

from nanopores.models import pughpore as pugh
from nanopores.tools.solvers import cache_forcefield
from nanopores.tools.utilities import collect_dict
from nanopores.models.diffusion import diffusivity

@cache_forcefield("diffz_pugh",
    default=dict(dim=2, h=1., Nmax=1e5, diamPore=6., rMolecule=0.11, H=100.))
def diff2D(X, **params):
    for x, result in collect_dict(X):
        params["x0"] = x
        setup = pugh.Setup(**params)
        D = diffusivity(setup)
        result.new = dict(D=D)
    return result

def diff_profile_z_pugh(a=-36., b=36., N=40, nproc=1, **params):
    X = [[0., 0., z] for z in np.linspace(a, b, N)]
    return diff2D(X, nproc=nproc, **params)

# REMARK: this could be the blueprint to a general "cache_functions" wrapper
def cache_pugh_diffusivity(geoname="pugh2", mode="coupled", **params):
    #name = "D%s-%s" % (geoname, mode)
    name = "D%s" %(geoname,)

    if not fields.exists(name, **params):
        if not "cheapest" in params:
            params["cheapest"] = True
        setup = pugh.Setup(x0=None, **params)
        r = setup.geop.rMolecule
        diamPore = setup.geop.diamPore

        if mode == "coupled":
            data_z = diff_profile_z_pugh(diamPore=diamPore)
            data_r = diff_profile_plane(r)
        elif mode == "simple":
            data_z = None
            data_r = diff_profile_plane(r)
        elif mode == "profile":
            data_z = diff_profile_z_pugh(diamPore=diamPore)
            data_r = diff_profile_trivial(r)
        else:
            raise NotImplementedError

        functions = diffusivity_field(setup, r, ddata_z=data_z, ddata_r=data_r,
                                      boundary="dnab", poreregion="poreregion")
        fields.save_functions(name, params, **functions)
        fields.update()

    return name

def cache_pugh_diffusivity_old(**params):
    "the function above applied to the pugh pore"
    if not fields.exists("Dpugh", **params):
        setup_params = dict(params)
        r = setup_params.pop("r")
        ddata_z = fields.get_fields("pugh_diff2D", rMolecule=r)
        ddata_r = fields.get_fields("pugh_diff3D", rMolecule=r, bulkbc=True)
        if not "cheapest" in setup_params:
            setup_params["cheapest"] = True
        setup = pugh.Setup(x0=None, **setup_params)
        functions = diffusivity_field(setup, r, ddata_r, ddata_z)
        fields.save_functions("Dpugh", params, **functions)
        fields.update()

def get_pugh_diffusivity(**params):
    cache_pugh_diffusivity(**params)
    functions, mesh = fields.get_functions("Dpugh2", **params)
    return functions, mesh

def cache_pugh_diffusivity_alt(**params):
    "the function above applied to the pugh pore"
    if not fields.exists("Dpugh_alt", **params):
        setup_params = dict(params)
        r = setup_params.pop("r")
        ddata_pore = fields.get_fields("pugh_diff_pore", rMolecule=r)
        ddata_bulk = fields.get_fields("pugh_diff_bulk", rMolecule=r)
        if not "cheapest" in setup_params:
            setup_params["cheapest"] = True
        setup = pugh.Setup(x0=None, **setup_params)
        functions = diffusivity_field_alt(setup, r, ddata_pore, ddata_bulk)
        fields.save_functions("Dpugh_alt", params, **functions)
        fields.update()

def get_pugh_diffusivity_alt(**params):
    cache_pugh_diffusivity_alt(**params)
    functions, mesh = fields.get_functions("Dpugh_alt", **params)
    return functions, mesh

if __name__ == "__main__":
    params = nanopores.user_params(dim=2, r=0.11, h=1., Nmax=1e5)
    functions = get_pugh_diffusivity(**params)
    D = functions["D"]
    D = dolfin.as_matrix(np.diag([D[i] for i in range(params.dim)]))
    dolfin.plot(functions["D"][0], title="Dx")
    dolfin.plot(functions["D"][1], title="Dz")
    dolfin.plot(functions["dist"], title="dist")
    dolfin.interactive()
