# (c) 2017 Gregor Mitscha-Baude
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
    dt = 10., # time step [ns]
    walldist = 1.6, # in multiples of radius, should be >= 1
    margtop = 50.,
    margbot = 20.,
    zstart = 70, # 46.5
    rstart = 70, # 42.
)

# domains are places where molecule can bind
# and/or be reflected after collision
#domain_params = dict(
#    cyl = False, # determines whether rz or xyz coordinates are passed to .inside
#    walldist = 1.5, # multiple of radius that determines what counts as collision
#
#    exclusion = True,
#    minsize = 0.01, # accuracy when performing reflection
#
#    binding = False,
#    eps = 1., # margin in addition to walldist, determines re-attempting
#    p = 0.1, # binding probability for one attempt
#    t = 1e6, # mean of exponentially distributed binding duration [ns]
#    dx = 0.4, # width of bond energy barrier [nm]
#    use_force = True, # if True, t_mean = t*exp(-|F|*dx/kT)
#)

pore = nanopores.get_pore(**params)
rw = randomwalk.RandomWalk(pore, **params)

receptor = randomwalk.Ball([39., 0., 46.5], 1.25)
rw.add_domain(receptor, exclusion=True, walldist=1.2,
              binding=True, eps=1., t=1e6, p=0.5)
#rw.add_wall_binding(t=1e4, p=0.1, eps=0.1)
randomwalk.run(rw, "rw_wei", a=-3, b=6, save_count=1000)