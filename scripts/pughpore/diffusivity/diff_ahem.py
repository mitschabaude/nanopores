# (c) 2017 Gregor Mitscha-Baude
import numpy as np
from math import sinh, acosh
from nanopores.models.diffusion_interpolation import diffusivity_field_simple
from nanopores.tools import fields

def Cp(h, r):
    x = r/h
    return 1. - 9./16.*x + 1./8.*x**3 - 45./256.*x**4 - 1/16.*x**5

def Cn(l, r):
    alpha = acosh(l/r)
    s = 0.
    for n in range(1, 100):
        n = float(n)
        K = n*(n+1)/(2*n-1)/(2*n+3)
        s += K*((2*sinh((2*n+1)*alpha)+(2*n+1)*sinh(2*alpha))/(4*(sinh((n+.5)*alpha))**2-(2*n+1)**2*(sinh(alpha))**2) - 1)
    return 1./((4./3.)*sinh(alpha)*s)

def matrix(d):
    return [[d[0], 0., 0.], [0., d[1], 0.], [0., 0., d[2]]]

def diff_profile_bulk(setup):
    geop = setup.geop

    r = geop.rMolecule
    eps = 1e-8
    x = np.linspace(r+eps, r*16., 100)

    Dn = np.array([Cn(xx, r) for xx in x])
    Dt = Cp(x, r)

    data = dict(x=list(x), D=map(matrix, zip(Dn, Dt, Dt)))
    return data

def simple_D(setup):
    r = setup.geop.rMolecule
    data = diff_profile_bulk(setup)
    functions = diffusivity_field_simple(setup, r, data)
    return functions["D"]

def cache_diffusivity(geoname="alphahem", **params):
    name = "D%s" % (geoname,)

    if not fields.exists(name, **params):
        setup = nanopore.Setup(**params)
        r = setup.geop.rMolecule
        data = diff_profile_bulk(setup)
        functions = diffusivity_field_simple(setup, r, data)
        fields.save_functions(name, params, **functions)
        fields.update()
        print

def get_diffusivity(geoname="alphahem", **params):
    cache_diffusivity(geoname, **params)
    name = "D%s" % (geoname,)
    functions, mesh = fields.get_functions(name=name, **params)
    D = functions["D"]
    return D

if __name__ == "__main__":
    import nanopores.models.nanopore as nanopore
    from matplotlib import pyplot as plt
    import dolfin, os
    fields.set_dir(os.path.expanduser("~") + "/Dropbox/nanopores/fields")

    params = dict(dim=2, Nmax=.1e5, h=1., ahemqsuniform=True, rMolecule=0.11)
    setup = nanopore.Setup(**params)
    data = diff_profile_bulk(setup)
    x = data["x"]
    Dt = [D[2][2] for D in data["D"]]
    Dn = [D[0][0] for D in data["D"]]

    plt.plot(x, Dt, ".-", label=r"$D_t$")
    plt.plot(x, Dn, ".-", label=r"$D_t$")
    plt.legend(loc="lower right")
#    plt.show()
#
#    D = get_diffusivity(**params)
#    dolfin.plot(D[0], title="Dx")
#    dolfin.plot(D[1], title="Dz")
#    dolfin.interactive()
    ddata = dict(name="Dalphahem", dim=2, Nmax=.4e5, h=1., ahemqsuniform=True, rMolecule=0.11)
    setup = nanopore.Setup(diffusivity_data=ddata, **params)
    _, pnps = nanopore.solve(setup, True)
    print nanopore.get_forces(setup, pnps)