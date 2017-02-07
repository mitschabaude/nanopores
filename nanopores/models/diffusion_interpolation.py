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
    D = np.diag([Dn, Dt, Dt]) if len(n)==3 else np.diag([Dn, Dt])
    U, S, V = np.linalg.svd(np.matrix(n))
    # U, S = 1., |n|
    # V = orthogonal coordinate basis with n/|n| as first vector
    D = np.dot(V, np.dot(D, V.T))
    return np.diag(np.diag(D))

def diffusivity_field(setup, r, ddata_r, ddata_z):
    "interpolates diffusivity field defined on the geometry given by setup"
    "needs: molecule radius r, radial and vertical D profiles ddata_r, ddata_z"
    "the geometry needs to have subdomains poreregion and bulkfluid"

    # build 1D interpolations from data
    data = ddata_z
    z = [x[2] for x in data["x"]]
    data, z = fields._sorted(data, z)
    Dz = data["D"]

    data = ddata_r
    x = [x[0] for x in data["x"]]
    data, x = fields._sorted(data, x)
    Dt = [d[2][2] for d in data["D"]]
    Dn = [d[0][0] for d in data["D"]]

    eps = 1e-2
    x = [-eps, r] + x + [100.]
    Dn = [d/Dn[-1] for d in Dn]
    Dt = [d/Dt[-1] for d in Dt]
    Dn = [0., eps] + Dn + [1.]
    Dt = [0., eps] + Dt + [1.]

    fz = interp1d(z, Dz)
    fn = interp1d(x, Dn)
    ft = interp1d(x, Dt)

    setup.prerefine(visualize=True)
    dist = distance_boundary_from_geo(setup.geo)
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

    D = lambda i: dict(
        bulkfluid = lambda x: DBulk(x, i),
        poreregion = lambda x: DPore(x, i),
    )

    dim = setup.geop.dim
    DD = [harmonic_interpolation(setup, subdomains=D(i)) for i in range(dim)]
    D = dolfin.Function(VV)
    dolfin.assign(D, DD)

    return dict(dist=dist, D=D)

from nanopores.models import pughpore as pugh

# REMARK: this could be the blueprint to a general "cache_functions" wrapper
def cache_pugh_diffusivity(**params):
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
    functions, mesh = fields.get_functions("Dpugh", **params)
    return functions

if __name__ == "__main__":
    params = nanopores.user_params(dim=2, r=0.11, h=1., Nmax=1e5)
    functions = get_pugh_diffusivity(**params)
    D = functions["D"]
    D = dolfin.as_matrix(np.diag([D[i] for i in range(params.dim)]))
    dolfin.plot(functions["D"][0], title="Dx")
    dolfin.plot(functions["D"][1], title="Dz")
    dolfin.plot(functions["dist"], title="dist")
    dolfin.interactive()
