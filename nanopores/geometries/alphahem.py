# (c) 2017 Gregor Mitscha-Baude
from alphahempoly import poly
from cylpore import get_geo as get_geo_generic

default = dict(
    dim = 2,
    geoname = "alphahem",
    porename = "alphahem",
    Htop = 7.,
    Hbot = 15.,
    R = 10.,
    cs=[-3, -6],
    zmem=-5.,
    #proteincs=[-5.],
)

def get_geo(h=1., recreate=False, **params):
    params = dict(default, **params)
    geo = get_geo_generic(poly, h=h, recreate=recreate, **params)
    return geo

if __name__ == "__main__":
    from nanopores import user_params
    params = user_params(default, x0=None, h=1.)
    geo = get_geo(recreate=True, **params)
    geo.plot_subdomains()
    geo.plot_boundaries(interactive=True)
