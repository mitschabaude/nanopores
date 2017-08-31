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
    dt = 10., # time step [ns]
    walldist = 1.5, # in multiples of radius, should be >= 1
    margtop = 50.,
    margbot = 20.,
    zstart = 60., # 46.5
    xstart = 38., # 42.
    rstart = 10.
)

#class RandomWalk(randomwalk.RandomWalk):
#    
#    # start at point [rstart, 0, zstart]
#    def initial(self):
#        rstart = self.params.rstart
#        zstart = self.params.zstart
#        if rstart is None:
#            rstart = 0.
#        if zstart is None:
#            zstart = self.ztop
#    
#        x = np.zeros((self.N, 3))
#        x[:, 0] = rstart
#        x[:, 1] = 0.
#        x[:, 2] = zstart
#        return x, np.abs(x[:, 0]), x[:, 2]

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

receptor = randomwalk.Ball([42. - 4., 0., 46.5], 1.25)
rw.add_domain(receptor, exclusion=True, walldist=1.2,
              binding=True, eps=1., t=1e6, p=0.5)
#rw.add_wall_binding(t=1e4, p=0.1, eps=0.1)
randomwalk.run(rw, "rw_wei", a=-3, b=1.5, save_count=1000)