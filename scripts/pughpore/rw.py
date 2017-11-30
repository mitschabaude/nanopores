# (c) 2017 Gregor Mitscha-Baude
"""Random walks in Pugh pore on 2D proxy geometry, to determine distribution of
attempt time needed to fit binding parameters."""
import numpy as np
import matplotlib.pyplot as plt
import nanopores
import nanopores.models.randomwalk as randomwalk
from nanopores.tools import fields
fields.set_dir_mega()

params = nanopores.user_params(
    # general params
    geoname = "wei",
    dim = 2,
    rMolecule = 1.25, # 6.
    h = 5.,
    Nmax = 1e5,
    Qmol = 2., #15.,
    bV = -0.2,
    dp = 26.,
    geop = dict(dp = 26.),
    posDTarget = True,

    # random walk params
    N = 100, # number of (simultaneous) random walks
    dt = .5, # time step [ns]
    walldist = 2., # in multiples of radius, should be >= 1
    margtop = 60.,
    margbot = 0.,
    #zstart = 46.5, # 46.5
    #xstart = 0., # 42.
    rstart = 30,
    initial = "sphere",

    # receptor params
    lbind = 4., # length of binding zone [nm]
)