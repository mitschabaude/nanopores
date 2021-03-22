# (c) 2017 Gregor Mitscha-Baude
"""
this module actually provides very general tools for diffusivity
interpolation as soon as models.nanopore.Setup can handle general geometries.
"""
import numpy as np
from nanopores.models.diffusion import diffusivity
from nanopores.models.diffusion_interpolation import (diffusivity_field,
    diff_profile_plane, diff_profile_trivial)
from nanopores.tools import fields, collect_dict
from nanopores.tools.solvers import cache_forcefield
from nanopores.models import nanopore

def simple_D(setup):
    r = setup.geop.rMolecule
    functions = diffusivity_field(setup, r, boundary="poresolidb")
    return functions["D"]

def cache_diffusivity_simple(geoname="alphahem", **params):
    name = "D%s" % (geoname,)

    if not fields.exists(name, **params):
        setup = nanopore.Setup(**params)
        r = setup.geop.rMolecule
        functions = diffusivity_field(setup, r, boundary="poresolidb")
        fields.save_functions(name, params, **functions)
        fields.update()

def get_diffusivity_simple(geoname="alphahem", **params):
    cache_diffusivity_simple(geoname, **params)
    name = "D%s" % (geoname,)
    functions, mesh = fields.get_functions(name=name, **params)
    D = functions["D"]
    return D

def cache_diffusivity(geoname="alphahem", mode="coupled", **params):
    name = "D%s-%s" % (geoname, mode)

    if not fields.exists(name, **params):
        setup = nanopore.Setup(**params)
        r = setup.geop.rMolecule

        if mode == "coupled":
            data_z = diff_profile_z_ahem(**params)
            data_r = diff_profile_plane(r)
        elif mode == "simple":
            data_z = None
            data_r = diff_profile_plane(r)
        elif mode == "profile":
            data_z = diff_profile_z_ahem(**params)
            data_r = diff_profile_trivial(r)
        else:
            raise NotImplementedError

        functions = diffusivity_field(setup, r, ddata_z=data_z, ddata_r=data_r,
                                      boundary="poresolidb", poreregion="pore")
        fields.save_functions(name, params, **functions)
        fields.update()

    return name

def get_diffusivity(geoname="alphahem", mode="coupled", **params):
    name = cache_diffusivity(geoname, mode, **params)
    functions, mesh = fields.get_functions(name=name, **params)
    D = functions["D"]
    return D

@cache_forcefield("diffz")
def diff2D(X, **params):
    for x, result in collect_dict(X):
        params["x0"] = x
        setup = nanopore.Setup(**params)
        D = diffusivity(setup)
        result.new = dict(D=D)
    return result

def diff_profile_z_ahem(a=-10.3, b=0.05, N=20, **params):
    X = [[0., 0., z] for z in np.linspace(a, b, N)]
    return diff2D(X, name="diffz-ahem", nproc=5, **params)

if __name__ == "__main__":
    from matplotlib import pyplot as plt
    import dolfin
    fields.set_dir_dropbox()

    #params = dict(dim=2, Nmax=.1e5, h=1., ahemqsuniform=True, rMolecule=0.11)
    params = dict(dim=2, Nmax=1e5, h=.5, ahemqsuniform=True, rMolecule=0.11)
    data = diff_profile_plane(r=params["rMolecule"])
    x = data["x"]
    Dt = [D[2][2] for D in data["D"]]
    Dn = [D[0][0] for D in data["D"]]

    plt.plot(x, Dt, ".-", label=r"$D_t$")
    plt.plot(x, Dn, ".-", label=r"$D_n$")
    plt.legend(loc="lower right")
    plt.show()
#
    data = diff_profile_z_ahem(**params)
    z = [x0[2] for x0 in data["x"]]
    Dz = data["D"]
    plt.figure()
    plt.plot(z, Dz, ".-")
    #plt.show()

    plotD = False
    solve = False

    if plotD:
        D = get_diffusivity(mode="profile", **params)
        dolfin.plot(D[0], title="Dx")
        dolfin.plot(D[1], title="Dz")
        dolfin.interactive()

    if solve:
        print("TRYING TO SOLVE")
        ddata = dict(params, name="Dalphahem")
        setup = nanopore.Setup(diffusivity_data=ddata, **params)
        _, pnps = nanopore.solve(setup, True)
        print(nanopore.get_forces(setup, pnps))