# (c) 2017 Gregor Mitscha-Baude
import numpy as np
import nanopores
import nanopores.models.randomwalk as randomwalk
from nanopores.tools import fields
fields.set_dir_dropbox()

params = nanopores.user_params(
    # general params
    geoname = "wei",
    dim = 2,
    rMolecule = 6.,
    h = 10.,
    Nmax = 1e5,
    Qmol = 15.,
    bV = -0.2,
    
    # random walk params
    N = 100, # number of (simultaneous) random walks
    dt = 1., # time step [ns]
    walldist = 1.5, # in multiples of radius, should be >= 1
    margtop = 50.,
    margbot = 20.,
    zstart = 60., # 46.5
    xstart = 38., # 42.
    rstart = 10.
)

receptor_params = dict(
    exclusion = False,
    walldist = 1.,
    #minsize = 0.01, # accuracy when performing reflection
    
    binding = True,
    eps = 0.1, # margin in addition to walldist, determines re-attempting [nm]
    t = 1e6, # mean of exponentially distributed binding duration [ns]
    p = 0.05, # binding probability for one attempt
    
    use_force = True, # if True, t_mean = t*exp(-|F|*dx/kT)
    dx = 0.1, # width of bond energy barrier [nm]
)

pore = nanopores.get_pore(**params)
rw = randomwalk.RandomWalk(pore, **params)

receptor = randomwalk.Ball([42. - 4., 0., 46.5], 1.25)
rw.add_domain(receptor, **receptor_params)
#rw.add_wall_binding(t=1e4, p=0.1, eps=0.1)
randomwalk.run(rw, "rw_wei", a=-3, b=1.5, save_count=1000)