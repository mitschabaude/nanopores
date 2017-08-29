# (c) 2017 Gregor Mitscha-Baude
import nanopores
from nanopores.models.randomwalk import (get_pore, RandomWalk, Ball, run)

params = nanopores.user_params(
    # general params
    geoname = "alphahem",
    dim = 2,
    rMolecule = .25,
    h = 1.,
    Nmax = 4e4,
    Qmol = 1.,
    bV = -0.4,
    # random walk params
    N = 100, # number of (simultaneous) random walks
    dt = 1e-2, # time step [ns]
    walldist = 2., # in multiples of radius, should be >= 1
    margtop = 2.,
    margbot = .5,
    rstart = 0.1,
    zstart = -5.,
)

pore = get_pore(**params)
rw = RandomWalk(pore, **params)
#receptor = Ball([1., 0., -30.], 8.)
#rw.add_domain(receptor, exclusion=True, walldist=1.,
#              binding=True, eps=1., t=1e6, p=0.1)
rw.add_wall_binding(t=1e4, p=0.1, eps=0.1)
run(rw, "rw_ahem", a=-3, b=6, save_count=1000)